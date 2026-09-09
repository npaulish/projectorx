""" Utility functions for ProjectorX. """

import os
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

import typer
from rich.console import Console

console = Console()

OPENMX_URL = "https://www.openmx-square.org/openmx3.9.tar.gz"
OPENMX_DIRNAME = "openmx3.9"

def ensure_openmx_exists(pao_path: Path, auto: bool = False):
    """
    Ensure OpenMX PAO directory exists.
    If not, optionally prompt the user or auto-download in CI.
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
            console.print(
                "[red]OpenMX is required. Please download manually from "
                f"{OPENMX_URL}[/red]"
            )
            raise typer.Exit(1)
    else:
        console.print(
            "[cyan]Running in non-interactive or CI mode — downloading automatically...[/cyan]"
            )


    # Use cached tarball if it exists
    if tar_path.exists():
        console.print(f"[blue]Using cached OpenMX tarball:[/blue] {tar_path}")
    else:
        console.print(f"[blue]Downloading OpenMX from {OPENMX_URL} ...[/blue]")
        try:
            urllib.request.urlretrieve(OPENMX_URL, tar_path)
        except Exception as e:
            console.print(f"[red]Failed to download OpenMX:[/red] {e}")
            raise typer.Exit(1) from e
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
        console.print(f"[red]Error:[/red] Expected PAO directory not found at {expected_pao}")
        raise typer.Exit(1)

    console.print("[green]✓ OpenMX setup complete![/green]")

    return expected_pao
