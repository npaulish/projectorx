"""JTH PAW unbound partial waves (``0S``/``0P``/``0D``, occupation 0) must not become projectors."""

import numpy as np
import pytest

from projectorx.upfdict import newUPFDict

R = np.exp(np.linspace(np.log(1e-3), np.log(15.0), 200))


def _chi(label: str, l: int, n: int, occupation: float, content: np.ndarray | None = None) -> dict:
    chi = {"l": l, "n": n, "occupation": occupation}
    chi["content"] = R * np.exp(-R) if content is None else content
    if label is not None:
        chi["label"] = label
    return chi


def _upf(chis: list[dict]) -> newUPFDict:
    return newUPFDict(
        "2.0.1",
        None,
        header={"element": "Zn", "has_so": False},
        mesh={"r": R},
        pswfc={"chi": chis},
    )


def _zn_like() -> newUPFDict:
    """Same labels and occupations as Zn.paw.pbesol.z_12.atompaw.jth.v1.1-std.upf."""
    unbound = np.sin(R)  # oscillating and non-decaying, like the real ones
    return _upf([
        _chi("4S", 0, 1, 2.0),
        _chi("0S", 0, 2, 0.0, unbound),
        _chi("0P", 1, 1, 0.0, unbound),
        _chi("0P", 1, 2, 0.0, unbound),
        _chi("3D", 2, 1, 10.0),
        _chi("0D", 2, 2, 0.0, unbound),
    ])


@pytest.mark.parametrize(
    ("chi", "expected"),
    [
        (_chi("0P", 1, 1, 0.0), True),  # JTH unbound partial wave
        (_chi("0P", 1, 1, "0.0000"), True),  # occupation given as text
        (_chi(" 0D", 2, 1, 0.0), True),  # padded label
        (_chi("3D", 2, 1, 10.0), False),  # bound, occupied
        (_chi("4P", 1, 1, 0.0), False),  # real but empty shell: keep it
        (_chi("0P", 1, 1, 2.0), False),  # 0-label but occupied: not ours to drop
        (_chi(None, 1, 1, 0.0), False),  # no label at all: never dropped
    ],
)
def test_is_unbound_partial_wave(chi, expected):
    assert newUPFDict.is_unbound_partial_wave(chi) is expected


def test_to_projectors_drops_unbound_only_when_asked():
    upf = _zn_like()
    assert [p.label for p in upf.to_projectors()] == ["4s", "0s", "0p", "0p", "3d", "0d"]
    assert [p.label for p in upf.to_projectors(drop_unbound=True)] == ["4s", "3d"]
    assert [c["label"] for c in upf.unbound_partial_waves()] == ["0S", "0P", "0P", "0D"]


def test_to_projectors_all_unbound_raises():
    upf = _upf([_chi("0S", 0, 1, 0.0), _chi("0P", 1, 1, 0.0)])
    with pytest.raises(ValueError, match="only unbound"):
        upf.to_projectors(drop_unbound=True)
