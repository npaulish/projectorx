# Usage

projectorx adds missing projectors to UPF pseudopotentials so they can be used
for Wannierization with projectability disentanglement (Y. Jiang et al.,
[arXiv:2507.06840](https://doi.org/10.48550/arXiv.2507.06840)). The result is a
`.dat` file in the format `pw2wannier90.x` reads for external projectors.

## Commands

### Extend one pseudopotential

```sh
projectorx extend-upf Fe.upf
```

| Option | Default | Meaning |
|---|---|---|
| `--output-dir DIR` | current directory | Where to write `<element>.dat`. Created if missing. |
| `--orbital 4p,5s` | (see below) | Orbitals to add. Repeat the option or separate with commas; case does not matter. |
| `--required-orbitals FILE` | built-in list | JSON file mapping each element to the orbitals it must have. |

### Extend a whole library

```sh
projectorx extend-library path/to/pseudos/ --output-dir external_projectors
```

Every `*.upf` / `*.UPF` file in the directory is processed. Files that fail are
skipped and listed at the end; the run continues. Outputs are named by element,
so if two files share an element, the first in sorted filename order is used
and the other is reported as skipped.

| Option | Default | Meaning |
|---|---|---|
| `--output-dir DIR` | `./external_projectors` | Where to write the outputs. Created if missing. |
| `--required-orbitals FILE` | built-in list | As above. |

## Which orbitals are added

Unless `--orbital` is given, projectorx adds every orbital that the
required-orbitals list names for the element and the pseudopotential lacks. The
built-in list ([`required_orbitals.json`](../src/projectorx/required_orbitals.json))
covers H to Rn, except the lanthanides (La to Lu); for example, Fe requires
`4s, 4p, 3d`. Elements not in the list (lanthanides, actinides, Fr, Ra) fail
unless you pass `--orbital` or your own `--required-orbitals` file.

Each missing orbital is built as a hydrogenic function. Only s, p and d
orbitals are supported, and its shape parameter alpha is chosen in one of two
ways:

- **The pseudopotential already has a shell with the same l** (e.g. it has 3p
  and 4p is missing): alpha is chosen so the new orbital is orthogonal to the
  shell directly below it (3p). If that shell is missing, or no orthogonal
  orbital exists for alpha in [1, 10), the file fails with an error.
- **The pseudopotential has no shell with that l at all** (e.g. K has no d,
  so 3d is missing): alpha is fitted to the matching pseudo-atomic orbital
  (PAO) from the OpenMX 3.9 library.

### OpenMX PAO files

The first time a PAO is needed, projectorx offers to download OpenMX 3.9
(<https://www.openmx-square.org>) and extracts it to
`~/.projectorx/paolibs/openmx3.9/`. Declining, or a failed download, stops the
run. If the environment variable `CI=true` is set, the download starts without
asking.

When OpenMX has several PAO files for an element, projectorx considers only the
standard ones (`<El><rc>.pao`, `<El><rc>S.pao`, `<El><rc>H.pao`) and picks:

1. the one whose `valence.electron` is closest to the UPF's `z_valence`, so both
   treat the same shells as valence;
2. on a tie, the one with the largest cutoff radius `rc`.

The chosen file is printed.

### JTH PAW datasets

JTH PAW files (made with atompaw) contain extra unbound partial waves labelled
`0S`, `0P`, `0D` with occupation 0. They are not orbitals, so projectorx leaves
them out of the output and prints which ones it dropped.

## Output

`<element>.dat` lists the radial functions on the pseudopotential's logarithmic
mesh, one column per projector:

```
<number of mesh points> <number of projectors>
<l of each projector>
<j of each projector>              # spin-orbit pseudopotentials only
<x> <r> <projector 1> <projector 2> ...
...
```

Spin-orbit pseudopotentials are supported: each new orbital gets one projector
per j = l ± 1/2.

`extend-library` also writes `projectors.json`, which lists for each element its
projectors' `label`, `l`, `j` (spin-orbit only) and `alpha`. `alpha` is `"UPF"`
for functions taken unchanged from the pseudopotential.

## Python API

Both commands can be called from Python:

```python
from pathlib import Path
from projectorx.cli import extend_upf_impl, extend_library_impl

dat_file = extend_upf_impl("Fe.upf", output_dir=Path("out"), orbitals=["4p"])

result = extend_library_impl(Path("pseudos"), output_dir=Path("external_projectors"))
result["projectors"]  # same content as projectors.json
result["failures"]  # {file name: error message}
```

If OpenMX is needed but unavailable, these raise
`projectorx.utils.OpenMXUnavailableError`.
