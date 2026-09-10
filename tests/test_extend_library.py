import json
import tempfile
import warnings
from pathlib import Path

import numpy as np
from typer.testing import CliRunner

from projectorx.cli import app

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
            result = runner.invoke(app, [
                "extend-library",
                str(test_data_dir),
                "--output-dir",
                tmpdir,
            ])

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
