"""Shell detection must ignore placeholder ``0<l>`` labels (JTH PAW unbound partial waves)."""

from unittest.mock import patch

import numpy as np
import pytest

from projectorx import utils
from projectorx.cli import _add_orbitals, _existing_shells
from projectorx.fit_hydrogenics import r_hydrogenic
from projectorx.projectors import newProjector, newProjectors

X = np.linspace(np.log(1e-4), np.log(40.0), 1500)  # logarithmic mesh, as in UPF files
R = np.exp(X)


def _proj(label: str, l: int, n: int, alpha: float) -> newProjector:
    """A hydrogenic stand-in for a pseudo-wavefunction (`n` = number of radial nodes)."""
    return newProjector(X, r_hydrogenic(R, l, n, alpha), l, label=label)


def _jth_zn_like() -> newProjectors:
    """Zn JTH PAW labels: 4S, 0S, 0P, 0P, 3D, 0D (no real p shell)."""
    return newProjectors(
        [
            _proj("4s", 0, 0, 2.0),
            _proj("0s", 0, 1, 2.0),
            _proj("0p", 1, 0, 2.0),
            _proj("0p", 1, 0, 3.0),
            _proj("3d", 2, 0, 3.0),
            _proj("0d", 2, 0, 5.0),
        ]
    )


def _jth_mo_like() -> newProjectors:
    """Mo JTH PAW labels: 4S, 5S, 0P, 4P, 4D, 0D (one real p shell, 4p)."""
    return newProjectors(
        [
            _proj("4s", 0, 0, 3.0),
            _proj("5s", 0, 1, 3.0),
            _proj("0p", 1, 0, 2.0),
            _proj("4p", 1, 0, 4.0),
            _proj("4d", 2, 0, 4.0),
            _proj("0d", 2, 0, 6.0),
        ]
    )


def test_existing_shells_ignores_zero_labels():
    assert _existing_shells(_jth_zn_like()) == ["3d", "4s"]
    assert _existing_shells(_jth_mo_like()) == ["4d", "4p", "4s", "5s"]


def test_missing_l_with_only_placeholders_uses_openmx():
    """Zn-like: 4p is missing and only ``0P`` placeholders exist, so OpenMX must be used.

    Before the fix this raised "Can't find inner projectors for 4p".
    """

    class ReachedOpenMX(Exception):
        pass

    with patch("projectorx.utils.ensure_openmx_exists", side_effect=ReachedOpenMX):
        with pytest.raises(ReachedOpenMX):
            _add_orbitals(_jth_zn_like(), "Zn", ["4p"], z_valence=12.0)


def test_placeholder_not_counted_as_inner_shell():
    """Mo-like: 5p is orthogonalised against 4p only.

    Before the fix ``0P`` + ``4P`` counted as two p shells, which asked for an
    unsupported hydrogenic function ("Only support n = 0 or 1 for p orbitals").
    """
    proj = _jth_mo_like()
    with patch.object(utils, "ensure_openmx_exists", side_effect=AssertionError("OpenMX not expected")):
        _add_orbitals(proj, "Mo", ["5p"], z_valence=14.0)

    added = proj[-1]
    assert (added.label, added.l) == ("5p", 1)
    inner_4p = next(p for p in proj if p.label == "4p")
    dr = np.diff(R, prepend=0.0)
    overlap = np.sum(added.y * inner_4p.y * dr)
    norm = np.sqrt(np.sum(added.y**2 * dr) * np.sum(inner_4p.y**2 * dr))
    assert abs(overlap / norm) < 0.05
