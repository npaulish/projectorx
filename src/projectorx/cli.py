"""Console script for projectorx."""

import os
from pathlib import Path
import re

import typer
from rich.console import Console


from projectorx import utils
from projectorx.fit_hydrogenics import fit_ortho_projectors, fit_rsq_projector, r_hydrogenic
from projectorx.projectors import newProjector, newProjectors
from projectorx.upfdict import newUPFDict

app = typer.Typer()
console = Console()

@app.command()
def main():
    """Console script for projectorx."""
    console.print("Replace this message by putting your code into "
               "projectorx.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")

@app.command("extend-upf")
def extend_upf(
    input_file: str = typer.Argument(..., help="Path to the input UPF file"),
    output_dir: str = typer.Argument(..., help="Path to the output directory"),
):
    """
    Extend a UPF pseudopotential with additional projectors.
    """
    console.print(f"[bold green]Extending UPF pseudopotential[/bold green]")
    console.print(f"Input file: {input_file}")
    console.print(f"Output directory: {output_dir}")

    upffile = Path(input_file)

    upfdict = newUPFDict.from_upf(upffile)
    try:
        element = upfdict["header"]["element"]
    except KeyError as exc:
        raise KeyError("Can not find element information from UPF file") from exc
    proj = upfdict.to_projectors()

    pswfcs = [_.label for _ in proj]
    # pswfcs = ["3s", "3p", "3p", "3d", "3d", "4s"]
    pswfcs = list(set(pswfcs))
    # pswfcs = ["3s", "3p", "3d", "4s"]

    ###-> Set additional orbitals we need <-###
    additional_orbitals = ["4p"]

    # Use fit_projector to get fine alpha, and add additional Projector
    str2l = {"s": 0, "p": 1, "d": 2, "f": 3}

    try:
        proj[0].j
    except AttributeError:
        spin_orbit = False
    else:
        spin_orbit = True
    # Add additional projectors
    for orb in additional_orbitals:
        print(orb)
        l = str2l[orb[1]]
        n = len([_ for _ in pswfcs if orb[1] in _])
        if n == 0:
            # no inner shell found, can only find orbitals from third-party PAO library
            ###-> Set OpenMX location <-###
            pao_path = utils.ensure_openmx_exists(Path("../../paolibs/openmx3.9/DFT_DATA19/PAO/"))
            pao_file = (
                pao_path
                / [
                    _
                    for _ in os.listdir(pao_path)
                    if re.match(rf"{element}[0-9]*\.0.*\.pao", _)
                ][0]
            )
            pao = newProjectors.from_pao(pao_file, n, l)[0]
            alpha = fit_rsq_projector(pao, n)
            print(pao_file, ":", element, n, l, alpha)
            x = proj[0].x
            r = proj[0].r
            y = r_hydrogenic(r, l, n, alpha)
            if spin_orbit:
                proj.add_projector_soc(newProjector(x, y, l, label=orb, alpha=alpha))
            else:
                proj.add_projector(newProjector(x, y, l, label=orb, alpha=alpha))
        else:
            if not spin_orbit:
                ref = None
                for p in proj:
                    if (int(p.label[0]) == int(orb[0]) - 1) and (
                        p.label[1].lower() == orb[1].lower()
                    ):
                        ref = p
            else:  # spin_orbit
                ref = []
                for p in proj:
                    if (int(p.label[0]) == int(orb[0]) - 1) and (
                        p.label[1].lower() == orb[1].lower()
                    ):
                        ref.append(p)

            if ref is None:
                raise ValueError(f"Cant find inner projectors for {orb}")
            if isinstance(ref, list):
                if not len(ref) in [1, 2]:
                    raise ValueError(
                        f"Wrong inner projectors for {orb}, found {len(ref)}"
                    )
            print("fit from pao ortho")
            if spin_orbit:
                for ref_ in ref:
                    print(n)
                    alpha = fit_ortho_projectors(ref_, n)
                    x = proj[0].x
                    r = proj[0].r
                    y = r_hydrogenic(r, l, n, alpha)
                    proj.add_projector(
                        newProjector(  # pylint: disable=unexpected-keyword-arg
                            x, y, l, j=ref_.j, label=orb, alpha=alpha
                        )
                    )
            else:
                alpha = fit_ortho_projectors(ref, n)
                x = proj[0].x
                r = proj[0].r
                y = r_hydrogenic(r, l, n, alpha)
                proj.add_projector(newProjector(x, y, l, label=orb, alpha=alpha))

    # Output .dat file
    proj.to_file(f"{output_dir}/{element}.dat")

if __name__ == "__main__":
    app()
