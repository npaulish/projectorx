"""Console script for projectorx."""

import os
from pathlib import Path
import re
import json

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
    console.print("[bold green]Extending UPF pseudopotential[/bold green]")
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
            console.print(f"[cyan]No inner shell found for {orb}, fitting from PAO file.[/cyan]")

            ###-> Set OpenMX location <-###
            pao_path = utils.ensure_openmx_exists(
                Path("~/.projectorx") / "paolibs/openmx3.9/DFT_DATA19/PAO/"
                )
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

@app.command("extend-family")
def extend_family(
    pseudo_family: str = typer.Argument(
        ...,
        help="AiiDA pseudo family to extend (label as set with aiida-pseudo package).",
    ),
    required_orbitals: Path = typer.Option(
        Path(__file__).parent / "required_orbitals.json",
        help="Path to a JSON file defining the required orbitals.",
    ),
    output_dir: Path = typer.Option(
        Path.cwd() / "external_projectors",
        help="Path to the output directory for results.",
    ),
):
    """
    Extend an existing AiiDA pseudo family with additional atomic projectors.
    """
    try:
        from aiida import load_profile, orm
    except ImportError as exc:
        console.print(
            "[red]Error:[/red] AiiDA is not installed. " \
            "This function is intended to provide additional orbitals for an AiiDA pseudo family" \
            " for Wannierisation with AiiDA workflows. " \
            "Please install with "
            "`pip install aiida-core`."
        )
        raise typer.Exit(1) from exc
    try:
        from aiida_wannier90_workflows.utils.pseudo.data import load_pseudo_metadata
    except ImportError as exc:
        console.print(
            "[red]Error:[/red] AiiDA-Wannier90-Workflows is not installed. Please install with "
            "`pip install aiida-wannier90-workflows`."
        )
        raise typer.Exit(1) from exc

    # Ensure OpenMX data
    pao_path = utils.ensure_openmx_exists(
        Path("~/.projectorx") / "paolibs/openmx3.9/DFT_DATA19/PAO/"
    )

    # Load AiiDA profile
    load_profile()

    # Load required orbitals and semicore data
    with open(required_orbitals, encoding="utf-8") as fp:
        required_orbital_list = json.load(fp)

    try:
        upfs = orm.load_group(pseudo_family)
    except KeyError as exc:
        console.print(
            f"[red]Error:[/red] Could not find pseudo family '{pseudo_family}'"
            " in the AiiDA database."
            " Try installing the family with `aiida-pseudo install` or check the family name."
        )
        raise typer.Exit(1) from exc

    try:
        pswfc = load_pseudo_metadata(f"semicore/{pseudo_family.replace('/', '_')}.json")
    except FileNotFoundError as exc:
        console.print(
            f"[red]Error:[/red] Could not find semicore files for pseudo family '{pseudo_family}'. "
            "This family is not supported by aiida-wannier90-workflows. "
            "Try using the latest version of aiida-wannier90-workflows package."
        )
        raise typer.Exit(1) from exc

    str2l = {"s": 0, "p": 1, "d": 2, "f": 3}
    projectors = {}

    output_dir = output_dir / pseudo_family.replace("/", "_")
    output_dir.mkdir(parents=True, exist_ok=True)

    for upf in upfs.nodes:
        element = upf.element
        if element not in required_orbital_list:
            continue

        console.print(f"[cyan]Processing element:[/cyan] {element}")

        pswfc[element]["additional"] = [
            ao
            for ao in required_orbital_list[element]
            if ao.upper() not in pswfc[element]["pswfcs"]
        ]

        upfdict = newUPFDict.from_str(upf.get_content())
        proj = upfdict.to_projectors()
        spin_orbit = hasattr(proj[0], "j")

        for addit_orb in pswfc[element]["additional"]:
            console.print(f"  → Adding orbital: {addit_orb}")
            l = str2l[addit_orb[1]]
            n = len([_ for _ in pswfc[element]["pswfcs"] if addit_orb[1].upper() in _])

            if n == 0:
                pao_file = next(
                    (pao_path / f for f in os.listdir(pao_path)
                     if re.match(rf"{element}[0-9]*\.0.*\.pao", f)),
                    None
                )
                if pao_file is None:
                    console.print(f"[red]No PAO file found for {element}![/red]")
                    continue

                pao = newProjectors.from_pao(pao_file, n, l)[0]
                alpha = fit_rsq_projector(pao, n)
                r, x = proj[0].r, proj[0].x
                y = r_hydrogenic(r, l, n, alpha)

                proj.add_projector_soc(newProjector(x, y, l, label=addit_orb, alpha=alpha)) if spin_orbit else proj.add_projector(
                    newProjector(x, y, l, label=addit_orb, alpha=alpha)
                )
            else:
                ref = [p for p in proj if (int(p.label[0]) == int(addit_orb[0]) - 1)
                       and (p.label[1].lower() == addit_orb[1].lower())]
                if not ref:
                    console.print(f"[red]No reference projector found for {addit_orb}![/red]")
                    continue

                alpha = fit_ortho_projectors(ref[0], n)
                r, x = proj[0].r, proj[0].x
                y = r_hydrogenic(r, l, n, alpha)
                proj.add_projector(newProjector(x, y, l, label=addit_orb, alpha=alpha))

        proj.to_file(output_dir / f"{element}.dat")

        projectors[element] = [
            {
                "label": p.label,
                "l": p.l,
                **({"j": p.j} if spin_orbit else {}),
                "alpha": p.alpha,
            }
            for p in proj
        ]

    with open(output_dir / "projectors.json", "w", encoding="utf-8") as fp:
        json.dump(projectors, fp, indent=2)

    console.print(f"[green]✓ Extended family saved to[/green] {output_dir}")

if __name__ == "__main__":
    app()
