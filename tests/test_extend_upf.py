from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from typer.testing import CliRunner

from projectorx.cli import _resolve_orbitals, app

runner = CliRunner()

DATA_DIR = Path(__file__).parent / "data"
SSSP_DIR = DATA_DIR / "sssp"
PAO_DIR = DATA_DIR / "pao"


def _assert_dat_matches(generated: Path, reference: Path, n_header: int) -> None:
    """Compare two projector ``.dat`` files: header lines exactly, data numerically."""
    gen_header = generated.read_text(encoding="utf-8").splitlines()[:n_header]
    ref_header = reference.read_text(encoding="utf-8").splitlines()[:n_header]
    assert gen_header == ref_header, f"Header mismatch: {gen_header} != {ref_header}"

    ref_data = np.loadtxt(reference, skiprows=n_header)
    gen_data = np.loadtxt(generated, skiprows=n_header)
    assert ref_data.shape == gen_data.shape, (
        f"Shape mismatch: reference {ref_data.shape}, generated {gen_data.shape}"
    )
    np.testing.assert_allclose(
        gen_data,
        ref_data,
        rtol=1e-5,
        atol=1e-7,
        err_msg=f"{generated.name} differs from {reference.name} beyond tolerance",
    )


@pytest.mark.parametrize(
    ("upf", "reference", "element", "n_header"),
    [
        # Spin-orbit: 4p orthogonalised against 3p (j=1/2 and j=3/2).
        (DATA_DIR / "Co.upf", DATA_DIR / "Co_reference.dat", "Co", 3),
        # No spin-orbit: 4p orthogonalised against 3p (SSSP, GBRV ultrasoft).
        (
            SSSP_DIR / "Ti.us.pbesol.z_12.uspp.gbrv.v1.4.upf",
            SSSP_DIR / "Ti_reference.dat",
            "Ti",
            2,
        ),
    ],
    ids=["Co-soc", "Ti"],
)
def test_extend_upf_orthogonalised(tmp_path, upf, reference, element, n_header):
    """`extend-upf` builds missing orbitals from an inner shell without touching OpenMX."""
    with patch(
        "projectorx.utils.ensure_openmx_exists",
        side_effect=AssertionError("OpenMX PAO path should not be used here"),
    ):
        result = runner.invoke(app, ["extend-upf", str(upf), "--output-dir", str(tmp_path)])

    assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception!r}"
    output_file = tmp_path / f"{element}.dat"
    assert output_file.exists(), "Output .dat file was not created"
    _assert_dat_matches(output_file, reference, n_header)


def test_extend_upf_creates_missing_output_dir(tmp_path):
    """A non-existent --output-dir (including parents) is created, as extend-library does."""
    output_dir = tmp_path / "new" / "nested"
    upf = SSSP_DIR / "Ti.us.pbesol.z_12.uspp.gbrv.v1.4.upf"

    result = runner.invoke(app, ["extend-upf", str(upf), "--output-dir", str(output_dir)])

    assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception!r}"
    _assert_dat_matches(output_dir / "Ti.dat", SSSP_DIR / "Ti_reference.dat", n_header=2)


def test_extend_upf_from_openmx_pao(tmp_path):
    """K has no d shell, so 3d must be fitted to the OpenMX PAO (K14.0.pao, 9 valence e-)."""
    upf = SSSP_DIR / "K.nc.pbesol.z_9.oncvpsp4.sg15.v0.upf"

    with patch("projectorx.utils.ensure_openmx_exists", return_value=PAO_DIR) as mock_ensure:
        result = runner.invoke(app, ["extend-upf", str(upf), "--output-dir", str(tmp_path)])

    assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception!r}"
    mock_ensure.assert_called_once()
    assert "K14.0.pao" in result.output.replace("\n", "")
    _assert_dat_matches(tmp_path / "K.dat", SSSP_DIR / "K_reference.dat", n_header=2)


def test_extend_upf_orbital_option_is_case_insensitive(tmp_path):
    """`--orbital 4P` used to fail with KeyError: 'P'."""
    upf = SSSP_DIR / "Ti.us.pbesol.z_12.uspp.gbrv.v1.4.upf"

    result = runner.invoke(app, ["extend-upf", str(upf), "--output-dir", str(tmp_path), "--orbital", "4P"])

    assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception!r}"
    _assert_dat_matches(tmp_path / "Ti.dat", SSSP_DIR / "Ti_reference.dat", n_header=2)


@pytest.mark.parametrize(
    ("given", "expected"),
    [
        (["4p"], ["4p"]),
        (["4P,5S"], ["4p", "5s"]),
        (["4p", "5s"], ["4p", "5s"]),
        ([" 4p , 3D ,"], ["4p", "3d"]),
    ],
)
def test_resolve_orbitals_override_normalised(given, expected):
    assert _resolve_orbitals(None, "Ti", given, {}) == expected


@pytest.mark.parametrize("bad", ["p4", "4x", "10p", "4", "4pp"])
def test_resolve_orbitals_override_invalid(bad):
    with pytest.raises(ValueError, match="Invalid orbital"):
        _resolve_orbitals(None, "Ti", [bad], {})
