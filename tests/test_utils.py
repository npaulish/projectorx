from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

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
    tar_path = tmp_path / "openmx3.9.tar.gz"
    tmp_tar_path = tmp_path / "openmx3.9.tar.gz.part"

    # mock urlretrieve to simulate downloading to the temporary ``.part`` file
    def fake_urlretrieve(url, filename):
        Path(filename).touch()
        return (str(filename), None)

    mock_urlretrieve.side_effect = fake_urlretrieve

    # mock tarfile to simulate extraction
    mock_tar = MagicMock()
    mock_tar_open.return_value.__enter__.return_value = mock_tar

    def fake_extractall(path, **kwargs):
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
    mock_urlretrieve.assert_called_once_with(utils.OPENMX_URL, tmp_tar_path)
    assert not tmp_tar_path.exists()  # renamed to tar_path, not left behind
    assert not tar_path.exists()  # cleaned up after extraction (non-CI mode)
    mock_tar.extractall.assert_called_once_with(tmp_path, filter="data")
    assert result.exists()
    assert "openmx3.9/DFT_DATA19/PAO" in str(result)


@patch("urllib.request.urlopen", side_effect=utils.urllib.error.URLError("offline"))
def test_ensure_openmx_exists_download_fails(mock_urlopen, tmp_path):
    """If download fails, the function should raise OpenMXUnavailableError."""
    pao_path = tmp_path / "openmx3.9/DFT_DATA19/PAO"

    with patch("urllib.request.urlretrieve", side_effect=Exception("Network error")):
        with pytest.raises(utils.OpenMXUnavailableError, match="Network error"):
            utils.ensure_openmx_exists(pao_path, auto=True)
    assert not (tmp_path / "openmx3.9.tar.gz.part").exists()


@patch("urllib.request.urlopen", side_effect=utils.urllib.error.URLError("offline"))
@patch("urllib.request.urlretrieve")
def test_ensure_openmx_exists_declined(mock_urlretrieve, mock_urlopen, tmp_path, monkeypatch):
    """Declining the download prompt raises OpenMXUnavailableError and downloads nothing."""
    monkeypatch.setenv("CI", "false")
    pao_path = tmp_path / "openmx3.9/DFT_DATA19/PAO"

    with patch("typer.confirm", return_value=False):
        with pytest.raises(utils.OpenMXUnavailableError, match="declined"):
            utils.ensure_openmx_exists(pao_path)
    mock_urlretrieve.assert_not_called()


def _write_pao(directory: Path, name: str, valence: float) -> None:
    (directory / name).write_text(f"AtomSpecies 26\nvalence.electron {valence}\n", encoding="utf-8")


@pytest.fixture(name="fe_pao_dir")
def fe_pao_dir_fixture(tmp_path):
    """A fake PAO directory mimicking the OpenMX Fe variants."""
    for name, valence in [
        ("Fe4.0H.pao", 16.0),
        ("Fe6.0H.pao", 16.0),
        ("Fe10.0H.pao", 16.0),
        ("Fe5.0S.pao", 14.0),
        ("Fe6.0S.pao", 14.0),
        ("Fe10.0S.pao", 14.0),
        ("Fe5.0p.pao", 14.0),  # non-standard variants must be ignored
        ("Fe12.0Hp.pao", 16.0),
        ("Fe12.0_OC.pao", 14.0),
        ("F7.0.pao", 7.0),  # different element sharing a prefix
    ]:
        _write_pao(tmp_path, name, valence)
    return tmp_path


@pytest.mark.parametrize(
    ("z_valence", "expected"),
    [
        (16.0, "Fe10.0H.pao"),  # exact match with hard variant, largest cutoff
        (14.0, "Fe10.0S.pao"),  # exact match with soft variant, largest cutoff
        (14.6, "Fe10.0S.pao"),  # closest valence wins
        (8.0, "Fe10.0S.pao"),  # no exact match: closest valence wins
        (None, "Fe10.0H.pao"),  # unknown valence: largest cutoff, then name
    ],
)
def test_select_pao_file(fe_pao_dir, z_valence, expected):
    """PAO selection is deterministic and prefers matching valence, then largest cutoff."""
    assert utils.select_pao_file(fe_pao_dir, "Fe", z_valence).name == expected


def test_select_pao_file_independent_of_listdir_order(fe_pao_dir):
    """The choice must not depend on the order os.listdir returns entries in."""
    names = sorted(p.name for p in fe_pao_dir.iterdir())
    picks = set()
    for order in (names, names[::-1]):
        with patch("projectorx.utils.os.listdir", return_value=order):
            picks.add(utils.select_pao_file(fe_pao_dir, "Fe", 14.0).name)
    assert picks == {"Fe10.0S.pao"}


def test_select_pao_file_missing_element(fe_pao_dir):
    """No standard file for the element raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError, match="Co"):
        utils.select_pao_file(fe_pao_dir, "Co", 17.0)
