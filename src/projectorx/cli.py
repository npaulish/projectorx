"""Console script for projectorx."""

import json
import os
import re
from pathlib import Path

import typer
from rich.console import Console

from projectorx import utils
from projectorx.fit_hydrogenics import fit_ortho_projectors, fit_rsq_projector, r_hydrogenic
from projectorx.projectors import newProjector, newProjectors
from projectorx.upfdict import newUPFDict

app = typer.Typer()
console = Console()

DEFAULT_REQUIRED_ORBITALS = Path(__file__).parent / "required_orbitals.json"


def _load_upf(upf_file: Path) -> newUPFDict:
    """Load a UPF file, tolerating the legacy ``Begin PSP_UPF`` / ``END_PSP`` wrapper some libraries use."""
    text = upf_file.read_text(encoding="utf-8")
    start = text.find("<UPF")
    end = text.rfind("</UPF>")
    if start != -1 and end != -1:
        text = text[start : end + len("</UPF>")]
    return newUPFDict.from_str(text)


def _get_element(upfdict: newUPFDict) -> str:
    """Read the element symbol from a UPF header.

    Works around two upf_tools/UPF-library quirks seen in real pseudopotential
    libraries: some files pad single-letter symbols with a space (e.g. ``"B "``),
    and upf_tools' attribute-type coercion turns a bare ``F`` (fluorine) into the
    Python boolean ``False``.
    """
    try:
        element = upfdict["header"]["element"]
    except KeyError as exc:
        raise KeyError("Can not find element information from UPF file") from exc
    if element is False:
        return "F"
    return str(element).strip()


def _resolve_orbitals(
    proj: newProjectors,
    element: str,
    orbitals_override: list[str] | None,
    required_orbital_list: dict,
) -> list[str]:
    """Work out which orbitals to add to `proj`.

    If `orbitals_override` is given, use it verbatim. Otherwise, fall back to
    whichever orbitals `required_orbital_list[element]` lists that aren't already
    present in `proj`.
    """
    if orbitals_override is not None:
        return [orb for group in orbitals_override for orb in group.split(",") if orb]

    pswfcs = list(set(p.label for p in proj))
    try:
        return [orb for orb in required_orbital_list[element] if orb not in pswfcs]
    except KeyError as exc:
        raise KeyError(
            f"No required orbitals defined for element '{element}'. "
            "Specify orbitals explicitly with --orbital."
        ) from exc


def _add_orbitals(proj: newProjectors, element: str, orbitals: list[str]) -> None:
    """Fit and add `orbitals` to `proj` in place."""
    pswfcs = list(set(p.label for p in proj))
    str2l = {"s": 0, "p": 1, "d": 2, "f": 3}
    spin_orbit = hasattr(proj[0], "j")

    for orb in orbitals:
        console.print(f"  → Adding orbital: {orb}")
        l = str2l[orb[1]]
        n = len([_ for _ in pswfcs if orb[1] in _])
        if n == 0:
            # no inner shell found, can only find orbitals from third-party PAO library
            console.print(f"[cyan]No inner shell found for {orb}, fitting from PAO file.[/cyan]")

            pao_path = utils.ensure_openmx_exists(
                Path.home() / ".projectorx" / "paolibs/openmx3.9/DFT_DATA19/PAO/"
            )
            pao_candidates = [
                f for f in os.listdir(pao_path) if re.match(rf"{element}[0-9]*\.0.*\.pao", f)
            ]
            if not pao_candidates:
                raise FileNotFoundError(f"No PAO file found for element {element} in {pao_path}")
            pao_file = pao_path / pao_candidates[0]
            pao = newProjectors.from_pao(pao_file, n, l)[0]
            alpha = fit_rsq_projector(pao, n)
            console.print(f"{pao_file} : {element} {n} {l} {alpha}")
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
                raise ValueError(f"Can't find inner projectors for {orb}")
            if isinstance(ref, list):
                if len(ref) not in [1, 2]:
                    raise ValueError(
                        f"Wrong inner projectors for {orb}, found {len(ref)}"
                    )
            if spin_orbit:
                for ref_ in ref:
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


@app.command("extend-upf")
def extend_upf(
    input_file: str = typer.Argument(..., help="Path to the input UPF file"),
    output_dir: Path = typer.Option(
        Path.cwd(),
        help="Path to the output directory for result.",
    ),
    required_orbitals: Path = typer.Option(
        DEFAULT_REQUIRED_ORBITALS,
        help="Path to a JSON file defining the required orbitals per element.",
    ),
    orbitals: list[str] = typer.Option(
        None,
        "--orbital",
        help=(
            "Orbital(s) to add, e.g. --orbital 4p,5s or --orbital 4p --orbital 5s. "
            "Defaults to whichever orbitals from --required-orbitals are missing from the UPF file."
        ),
    ),
):
    """
    Extend a UPF pseudopotential with additional projectors.
    """
    console.print("[bold green]Extending UPF pseudopotential[/bold green]")
    console.print(f"Input file: {input_file}")
    console.print(f"Output directory: {output_dir}")

    upfdict = _load_upf(Path(input_file))
    element = _get_element(upfdict)
    proj = upfdict.to_projectors()

    with open(required_orbitals, encoding="utf-8") as fp:
        required_orbital_list = json.load(fp)
    orbitals_to_add = _resolve_orbitals(proj, element, orbitals, required_orbital_list)

    if not orbitals_to_add:
        console.print(f"[green]✓ {element} already has all required orbitals, nothing to add.[/green]")
    _add_orbitals(proj, element, orbitals_to_add)

    # Output .dat file
    proj.to_file(f"{output_dir}/{element}.dat")


@app.command("extend-library")
def extend_library(
    input_dir: Path = typer.Argument(
        ...,
        help="Directory containing UPF pseudopotential files (*.upf / *.UPF).",
    ),
    output_dir: Path = typer.Option(
        Path.cwd() / "external_projectors",
        help="Path to the output directory for results.",
    ),
    required_orbitals: Path = typer.Option(
        DEFAULT_REQUIRED_ORBITALS,
        help="Path to a JSON file defining the required orbitals per element.",
    ),
):
    """
    Extend every UPF pseudopotential in a directory with additional projectors.

    For each file, adds whichever orbitals --required-orbitals lists for its
    element that aren't already present.
    """
    upf_files = sorted({*input_dir.glob("*.upf"), *input_dir.glob("*.UPF")})
    if not upf_files:
        console.print(f"[red]No .upf files found in {input_dir}[/red]")
        raise typer.Exit(1)

    with open(required_orbitals, encoding="utf-8") as fp:
        required_orbital_list = json.load(fp)

    output_dir.mkdir(parents=True, exist_ok=True)

    projectors_summary = {}
    failures = {}

    for upf_file in upf_files:
        console.print(f"[bold cyan]Processing:[/bold cyan] {upf_file.name}")
        try:
            upfdict = _load_upf(upf_file)
            element = _get_element(upfdict)
            proj = upfdict.to_projectors()

            orbitals_to_add = _resolve_orbitals(proj, element, None, required_orbital_list)
            if not orbitals_to_add:
                console.print(f"[green]✓ {element} already has all required orbitals.[/green]")
            _add_orbitals(proj, element, orbitals_to_add)
        except typer.Abort:
            # The user declined (or couldn't answer) the OpenMX download prompt;
            # every later element needing PAO fitting would hit the same wall.
            console.print("[red]Aborted: OpenMX download declined or unavailable.[/red]")
            raise
        except Exception as exc:  # noqa: BLE001 -- keep processing the rest of the library
            message = str(exc) or type(exc).__name__
            console.print(f"[red]✗ Failed to extend {upf_file.name}: {message}[/red]")
            failures[upf_file.name] = message
            continue

        proj.to_file(output_dir / f"{element}.dat")
        spin_orbit = hasattr(proj[0], "j")
        projectors_summary[element] = [
            {
                "label": p.label,
                "l": p.l,
                **({"j": p.j} if spin_orbit else {}),
                "alpha": p.alpha,
            }
            for p in proj
        ]

    with open(output_dir / "projectors.json", "w", encoding="utf-8") as fp:
        json.dump(projectors_summary, fp, indent=2)

    console.print(
        f"[green]✓ Extended {len(projectors_summary)}/{len(upf_files)} pseudopotentials, "
        f"saved to {output_dir}[/green]"
    )
    if failures:
        console.print(f"[yellow]⚠ {len(failures)} pseudopotential(s) failed:[/yellow]")
        for name, err in failures.items():
            console.print(f"    {name}: {err}")

@app.command("extend-family")
def extend_family(
    pseudo_family: str = typer.Argument(
        ...,
        help="AiiDA pseudo family to extend (label as set with aiida-pseudo package).",
    ),
    required_orbitals: Path = typer.Option(
        DEFAULT_REQUIRED_ORBITALS,
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
        Path.home() / ".projectorx" / "paolibs/openmx3.9/DFT_DATA19/PAO/"
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

                new_proj = newProjector(x, y, l, label=addit_orb, alpha=alpha)
                if spin_orbit:
                    proj.add_projector_soc(new_proj)
                else:
                    proj.add_projector(new_proj)
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
