import numpy as np
from pathlib import Path
import tempfile

import pytest
from typer.testing import CliRunner

from projectorx.cli import app

runner = CliRunner()

@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"

def test_extend_upf_creates_expected_output(test_data_dir):
    """Test that `extend-upf` produces the correct .dat output for Co.upf."""
    co_upf = test_data_dir / "Co.upf"
    reference = test_data_dir / "Co_reference.dat"

    # Create temporary directory for output
    with tempfile.TemporaryDirectory() as tmpdir:
        result = runner.invoke(app, [
            "extend-upf",
            str(co_upf),
            tmpdir,
        ])

        # Check CLI exited successfully
        assert result.exit_code == 0, f"CLI failed: {result.output}"

        # Check output file exists
        output_file = Path(tmpdir) / "Co.dat"
        assert output_file.exists(), "Output .dat file was not created"

        # Compare to reference file
        ref_data = np.loadtxt(reference, skiprows=3)
        gen_data = np.loadtxt(output_file, skiprows=3)

        # Check shapes
        assert ref_data.shape == gen_data.shape, (
            f"Shape mismatch: reference {ref_data.shape}, generated {gen_data.shape}"
        )

        # Numerically tolerant comparison
        np.testing.assert_allclose(
            gen_data,
            ref_data,
            rtol=1e-5,  # relative tolerance (5 parts in 1e5)
            atol=1e-7,  # absolute tolerance
            err_msg="Generated Co.dat differs from reference Co_reference.dat beyond tolerance",
        )
