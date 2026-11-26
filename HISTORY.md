# History

## Unreleased

Changes to results:

* Unbound partial waves in JTH PAW datasets (`0S`, `0P`, `0D` with occupation 0)
  are no longer written as projectors. This changes the output for elements with
  JTH pseudopotentials (e.g. As, B, Cd, Ga, Mo, Nb, Sb, Sn, Xe, Zn).
* The OpenMX PAO file used for fitting is now chosen deterministically (closest
  `valence.electron` to the UPF's `z_valence`, then largest cutoff). Before, it
  depended on directory listing order, so results for elements fitted to OpenMX
  (e.g. K, Cs, Ba) could differ between machines.

Fixes:

* Zn and Mo JTH PAW pseudopotentials no longer fail: `0`-labelled functions are
  not counted as existing shells.
* Failed or non-physical fits now raise an error instead of writing an all-zero
  projector, or silently using the edge of the alpha range.
* `extend-library` no longer lets two files for the same element overwrite each
  other: the first in sorted order is used and the rest are reported.
* Declining or failing the OpenMX download stops `extend-library` once, instead
  of prompting again for every element.
* `extend-upf` creates `--output-dir` if it does not exist.
* `--orbital` is case-insensitive and rejects malformed values.
* Spin-orbit `j` must be l ± 1/2.
* OpenMX is extracted with tarfile's `data` filter where available.

Other:

* `extend_upf_impl` and `extend_library_impl` raise
  `projectorx.utils.OpenMXUnavailableError` instead of `typer.Exit`.
* More tests (K, Ti, JTH-style labels, fitting, spin-orbit `j`); tests are no
  longer included in the source distribution.

## 0.1.0 (2025-10-08)

* First release on PyPI.
