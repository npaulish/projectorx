import json
import shutil
import tempfile
import warnings
from pathlib import Path
from unittest.mock import patch

import numpy as np
from typer.testing import CliRunner

from projectorx import utils
from projectorx.cli import app, extend_library_impl

runner = CliRunner()


def test_extend_library_processes_all_and_skips_broken_files():
    """extend-library should extend every valid UPF in a directory and skip broken ones without aborting.

    Reuses tests/data/Co.upf (already used by test_extend_upf.py) rather than
    shipping a whole pseudopotential library just for this test, plus a tiny
    garbage file to exercise the per-file failure handling.
    """
    test_data_dir = Path(__file__).parent / "data"
    reference = test_data_dir / "Co_reference.dat"

    with tempfile.TemporaryDirectory() as tmpdir:
        with warnings.catch_warnings():
            # upf_tools warns when it can't determine broken.upf's version.
            warnings.simplefilter("ignore")
            result = runner.invoke(
                app,
                [
                    "extend-library",
                    str(test_data_dir),
                    "--output-dir",
                    tmpdir,
                ],
            )

        assert result.exit_code == 0, f"CLI failed: {result.output}"

        # The one valid pseudopotential (Co.upf) should be extended correctly.
        output_file = Path(tmpdir) / "Co.dat"
        assert output_file.exists(), "Output .dat file for Co was not created"

        ref_data = np.loadtxt(reference, skiprows=3)
        gen_data = np.loadtxt(output_file, skiprows=3)
        assert ref_data.shape == gen_data.shape
        np.testing.assert_allclose(gen_data, ref_data, rtol=1e-5, atol=1e-7)

        # broken.upf should be recorded as a failure, not silently dropped or fatal.
        summary = json.loads((Path(tmpdir) / "projectors.json").read_text())
        assert list(summary.keys()) == ["Co"]


def test_extend_library_stops_once_when_openmx_unavailable(tmp_path):
    """If OpenMX is declined/unavailable, extend-library must stop instead of re-prompting per element."""
    data_dir = Path(__file__).parent / "data"
    input_dir = tmp_path / "pseudos"
    input_dir.mkdir()
    shutil.copy(data_dir / "Co.upf", input_dir / "Co.upf")
    shutil.copy(data_dir / "Co.upf", input_dir / "Co_copy.upf")

    # Co.upf has no f shell, so requiring 4f forces the OpenMX PAO path.
    required = tmp_path / "required.json"
    required.write_text(json.dumps({"Co": ["4f"]}), encoding="utf-8")

    with patch(
        "projectorx.utils.ensure_openmx_exists",
        side_effect=utils.OpenMXUnavailableError("OpenMX download declined."),
    ) as mock_ensure:
        result = runner.invoke(
            app,
            [
                "extend-library",
                str(input_dir),
                "--output-dir",
                str(tmp_path / "out"),
                "--required-orbitals",
                str(required),
            ],
        )

    assert result.exit_code == 1, result.output
    mock_ensure.assert_called_once()
    assert "declined" in result.output


def test_extend_library_duplicate_element_is_not_overwritten(tmp_path):
    """Two files for the same element: the first (sorted) is used, the second is reported.

    Outputs are named <element>.dat, so previously the second file silently
    overwrote the first one's output and summary entry.
    """
    data_dir = Path(__file__).parent / "data"
    input_dir = tmp_path / "pseudos"
    input_dir.mkdir()
    shutil.copy(data_dir / "Co.upf", input_dir / "Co.upf")
    shutil.copy(data_dir / "Co.upf", input_dir / "Co_other.upf")
    output_dir = tmp_path / "out"

    with patch(
        "projectorx.utils.ensure_openmx_exists",
        side_effect=AssertionError("OpenMX not expected"),
    ):
        result = extend_library_impl(input_dir, output_dir)

    assert list(result["projectors"]) == ["Co"]
    assert list(result["failures"]) == ["Co_other.upf"]
    assert "duplicate element Co" in result["failures"]["Co_other.upf"]
    assert "Co.upf" in result["failures"]["Co_other.upf"]

    ref = np.loadtxt(data_dir / "Co_reference.dat", skiprows=3)
    np.testing.assert_allclose(np.loadtxt(output_dir / "Co.dat", skiprows=3), ref, rtol=1e-5, atol=1e-7)
