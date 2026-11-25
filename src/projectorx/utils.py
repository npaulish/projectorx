""" Utility functions for ProjectorX. """

import os
import re
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

import typer
from rich.console import Console

console = Console()

OPENMX_URL = "https://www.openmx-square.org/openmx3.9.tar.gz"
OPENMX_DIRNAME = "openmx3.9"


class OpenMXUnavailableError(RuntimeError):
    """The OpenMX PAO library is missing and could not (or may not) be downloaded."""


def ensure_openmx_exists(pao_path: Path, auto: bool = False):
    """
    Ensure OpenMX PAO directory exists.
    If not, optionally prompt the user or auto-download in CI.

    :raises OpenMXUnavailableError: the user declined the download, the download
        failed, or the extracted archive does not contain the PAO directory.
    """
    if pao_path.exists():
        console.print(f"[green]✓ OpenMX PAO directory found:[/green] {pao_path}")
        return pao_path

    console.print(f"[yellow]⚠ OpenMX not found at {pao_path}[/yellow]")

    dest_dir = pao_path.parent.parent.parent  # e.g., .../openmx3.9/DFT_DATA19/PAO/
    tar_path = dest_dir / "openmx3.9.tar.gz"
    tar_path.parent.mkdir(parents=True, exist_ok=True)

    size_mb = None
    try:
        with urllib.request.urlopen(OPENMX_URL) as response:
            file_size = response.getheader("Content-Length")
            if file_size:
                size_mb = int(file_size) / (1024 * 1024)
    except urllib.error.URLError:
        console.print("[orange]Failed to connect to OpenMX URL.[/orange]")

    if size_mb:
        console.print(f"[blue]Approximate download size:[/blue] {size_mb:.1f} MB")
    else:
        console.print("[yellow]Approximate download size could not be determined.[/yellow]")

    ci_mode = os.getenv("CI", "false").lower() == "true"

    if not auto and not ci_mode:
        if not typer.confirm("Do you want to download and extract OpenMX automatically?"):
            raise OpenMXUnavailableError(
                f"OpenMX download declined. Download it manually from {OPENMX_URL} "
                f"and extract it so that {pao_path} exists."
            )
    else:
        console.print(
            "[cyan]Running in non-interactive or CI mode — downloading automatically...[/cyan]"
            )


    # Use cached tarball if it exists
    if tar_path.exists():
        console.print(f"[blue]Using cached OpenMX tarball:[/blue] {tar_path}")
    else:
        console.print(f"[blue]Downloading OpenMX from {OPENMX_URL} ...[/blue]")
        # Download to a temporary name and rename on success, so an interrupted
        # download is never mistaken for a complete, valid cached tarball later.
        tmp_path = tar_path.with_name(tar_path.name + ".part")
        try:
            urllib.request.urlretrieve(OPENMX_URL, tmp_path)
        except Exception as e:
            tmp_path.unlink(missing_ok=True)
            raise OpenMXUnavailableError(f"Failed to download OpenMX from {OPENMX_URL}: {e}") from e
        tmp_path.rename(tar_path)
        console.print(f"[green]✓ Downloaded[/green] {tar_path}")

    console.print("[blue]Extracting...[/blue]")

    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(dest_dir)

    # Only delete tarball if not in CI mode (to allow caching)
    if not ci_mode:
        tar_path.unlink(missing_ok=True)
        console.print("[blue]Cleaned up the tarball[/blue]")
    else:
        console.print(f"[blue]Preserving the tarball for caching:[/blue] {tar_path}")

    console.print(f"[green]✓ Extracted OpenMX into[/green] {dest_dir / OPENMX_DIRNAME}")

    expected_pao = dest_dir / OPENMX_DIRNAME / "DFT_DATA19" / "PAO"
    if not expected_pao.exists():
        raise OpenMXUnavailableError(f"Expected PAO directory not found at {expected_pao}")

    console.print("[green]✓ OpenMX setup complete![/green]")

    return expected_pao


_STANDARD_PAO_RE = r"{element}(?P<cutoff>\d+(?:\.\d+)?)(?P<variant>[SH]?)\.pao"


def _read_pao_valence(pao_file: Path) -> float | None:
    """Return the ``valence.electron`` count from an OpenMX PAO file header, if present."""
    with open(pao_file, encoding="utf-8", errors="replace") as fp:
        for line in fp:
            if line.strip().startswith("valence.electron"):
                try:
                    return float(line.split()[1])
                except (IndexError, ValueError):
                    return None
    return None


def select_pao_file(pao_dir: Path, element: str, z_valence: float | None = None) -> Path:
    """Deterministically pick the OpenMX PAO file to fit `element`'s missing orbitals to.

    Only the standard variants ``<El><rc>.pao``, ``<El><rc>S.pao`` (soft) and
    ``<El><rc>H.pao`` (hard) are considered; specialised ones (``p``, ``Hp``,
    ``_3d``, ``_OC``, ...) are ignored. Candidates are ranked by:

    1. how close the PAO's ``valence.electron`` is to the UPF's `z_valence`, so the
       PAO treats the same shells as valence and its first orbital of a given l is
       the same shell as the one missing from the UPF;
    2. largest cutoff radius ``rc`` (most complete radial tail);
    3. file name, as a final tie-break.

    :raises FileNotFoundError: no standard PAO file exists for `element`.
    """
    pattern = re.compile(_STANDARD_PAO_RE.format(element=re.escape(element)))
    candidates = []
    for name in os.listdir(pao_dir):
        match = pattern.fullmatch(name)
        if not match:
            continue
        path = Path(pao_dir) / name
        valence = _read_pao_valence(path)
        if z_valence is None:
            valence_mismatch = 0.0
        elif valence is None:
            valence_mismatch = float("inf")
        else:
            valence_mismatch = abs(valence - z_valence)
        candidates.append((valence_mismatch, -float(match["cutoff"]), name, path))

    if not candidates:
        raise FileNotFoundError(f"No PAO file found for element {element} in {pao_dir}")
    return min(candidates)[-1]
