from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import typer

from projectorx import utils  # adjust to your actual import path


def test_ensure_openmx_exists_already_present(tmp_path):
    """If PAO path exists, it should return immediately without downloading."""
    pao_path = tmp_path / "openmx3.9/DFT_DATA19/PAO"
    pao_path.mkdir(parents=True)

    result = utils.ensure_openmx_exists(pao_path)

    assert result == pao_path
    # Nothing should be downloaded or extracted


@patch("urllib.request.urlretrieve")
@patch("tarfile.open")
def test_ensure_openmx_exists_downloads_and_extracts(mock_tar_open, mock_urlretrieve, tmp_path, monkeypatch):
    """Simulate missing OpenMX directory and ensure download + extract workflow."""
    # Pin non-CI behavior (tarball cleanup) regardless of the ambient environment
    # -- GitHub Actions sets CI=true, which intentionally preserves the tarball.
    monkeypatch.setenv("CI", "false")

    pao_path = tmp_path / "openmx3.9/DFT_DATA19/PAO"

    # mock urlretrieve to simulate successful download
    mock_urlretrieve.return_value = (str(tmp_path / "openmx3.9.tar.gz"), None)

    # mock tarfile to simulate extraction
    mock_tar = MagicMock()
    mock_tar_open.return_value.__enter__.return_value = mock_tar

    def fake_extractall(path):
        # Simulate extraction by creating the expected PAO directory
        extracted_pao = Path(path) / "openmx3.9/DFT_DATA19/PAO"
        extracted_pao.mkdir(parents=True, exist_ok=True)

    mock_tar.extractall.side_effect = fake_extractall

    with patch("urllib.request.urlopen") as mock_urlopen:
        mock_resp = MagicMock()
        mock_resp.getheader.return_value = "12345678"
        mock_urlopen.return_value.__enter__.return_value = mock_resp

        result = utils.ensure_openmx_exists(pao_path, auto=True)

    # Assertions
    mock_urlretrieve.assert_called_once_with(utils.OPENMX_URL, tmp_path / "openmx3.9.tar.gz")
    mock_tar.extractall.assert_called_once_with(tmp_path)
    assert result.exists()
    assert "openmx3.9/DFT_DATA19/PAO" in str(result)


def test_ensure_openmx_exists_download_fails(tmp_path):
    """If download fails, the function should raise typer.Exit."""
    pao_path = tmp_path / "openmx3.9/DFT_DATA19/PAO"

    with patch("urllib.request.urlretrieve", side_effect=Exception("Network error")):
        with pytest.raises(typer.Exit):
            utils.ensure_openmx_exists(pao_path, auto=True)
