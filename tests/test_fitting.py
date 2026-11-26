"""Fits must fail loudly instead of returning a meaningless alpha, and j must be l +/- 1/2."""

from unittest.mock import patch

import numpy as np
import pytest

from projectorx.fit_hydrogenics import fit_ortho_projectors, fit_rsq_projector, hydrogenic, r_hydrogenic
from projectorx.projectors import newProjector

X = np.linspace(np.log(1e-4), np.log(40.0), 1500)  # logarithmic mesh, as in UPF files
R = np.exp(X)


def _hydrogenic(l: int, n: int, alpha: float) -> newProjector:
    return newProjector(X, r_hydrogenic(R, l, n, alpha), l)


# --- fit_rsq_projector (OpenMX PAO path) -------------------------------------------


@pytest.mark.parametrize(("l", "alpha"), [(0, 2.0), (1, 3.5), (2, 2.6)])
@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_fit_rsq_recovers_alpha(l, alpha, sign):
    """A hydrogenic target (of either sign) gives back its own alpha.

    OpenMX PAO files store the radial function R(r) itself (|R| ~ r^l near 0),
    not r * R(r) as UPF files do, so the target is built the same way.
    """
    target = newProjector(X, sign * hydrogenic(R, l, 0, alpha), l)
    assert fit_rsq_projector(target, 0) == pytest.approx(alpha, rel=1e-4)


def test_fit_rsq_both_fits_fail_raises():
    """Previously returned alpha = 0, i.e. an all-zero projector."""
    with patch("projectorx.fit_hydrogenics.curve_fit", side_effect=RuntimeError("no convergence")):
        with pytest.raises(RuntimeError, match="no convergence"):
            fit_rsq_projector(_hydrogenic(1, 0, 3.0), 0)


def test_fit_rsq_non_physical_alpha_raises():
    with patch(
        "projectorx.fit_hydrogenics.curve_fit",
        return_value=(np.array([-2.0]), np.array([[0.1]])),
    ):
        with pytest.raises(RuntimeError, match="non-physical alpha"):
            fit_rsq_projector(_hydrogenic(1, 0, 3.0), 0)


def test_fit_rsq_one_failed_sign_is_ignored():
    """Previously a failed sign got a fake variance of 100 and alpha = 0, which won
    whenever the successful fit's variance was larger than 100."""
    results = iter([RuntimeError("no convergence"), (np.array([3.0]), np.array([[1e3]]))])

    def fake_curve_fit(*args, **kwargs):
        result = next(results)
        if isinstance(result, Exception):
            raise result
        return result

    with patch("projectorx.fit_hydrogenics.curve_fit", side_effect=fake_curve_fit):
        assert fit_rsq_projector(_hydrogenic(1, 0, 3.0), 0) == 3.0


# --- fit_ortho_projectors (orthogonalisation path) ---------------------------------


@pytest.mark.parametrize(("l", "alpha_inner"), [(0, 3.0), (1, 4.0)])
def test_fit_ortho_gives_orthogonal_orbital(l, alpha_inner):
    inner = _hydrogenic(l, 0, alpha_inner)
    alpha = fit_ortho_projectors(inner, 1)
    outer = r_hydrogenic(R, l, 1, alpha)
    dr = np.diff(R, prepend=0.0)
    overlap = np.sum(inner.y * outer * dr) / np.sqrt(np.sum(inner.y**2 * dr) * np.sum(outer**2 * dr))
    assert abs(overlap) < 0.05


def test_fit_ortho_without_sign_change_raises():
    """A very compact inner p function: the overlap is positive for every alpha in
    [1, 10), so no orthogonal orbital exists there. Previously the grid edge was
    returned silently."""
    compact = newProjector(X, R**2 * np.exp(-50.0 * R), 1)
    with pytest.raises(ValueError, match="never changes sign"):
        fit_ortho_projectors(compact, 1)


def test_fit_ortho_zero_inner_raises():
    with pytest.raises(ValueError, match="is zero"):
        fit_ortho_projectors(newProjector(X, np.zeros_like(R), 1), 1)


# --- newProjector.j -----------------------------------------------------------------


@pytest.mark.parametrize(("l", "j"), [(0, 0.5), (1, 0.5), (1, 1.5), (2, 1.5), (2, 2.5), (3, 3.5)])
def test_valid_j(l, j):
    assert newProjector(X, R, l, j=j).j == j


@pytest.mark.parametrize(("l", "j"), [(1, 1.0), (0, -0.5), (2, 0.5), (1, 2.5), (2, 2.0)])
def test_invalid_j_raises(l, j):
    """j = l was previously accepted (any |l - j| < 1/2 passed)."""
    with pytest.raises(ValueError, match=f"Invalid j={j} for l={l}"):
        newProjector(X, R, l, j=j)


def test_no_j_means_no_spin_orbit():
    assert not hasattr(newProjector(X, R, 1), "j")
