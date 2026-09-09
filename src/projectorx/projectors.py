"""Override of upf_tools.Projector and .Projectors."""

from dataclasses import InitVar, dataclass, field
from pathlib import Path

import numpy as np
from upf_tools.projectors import Projector, Projectors

num2label = {0: "s", 1: "p", 2: "d", 3: "f"}


@dataclass
class newProjector(Projector):
    """Override of upf_tools.Projctor."""

    _x_min: float = field(default=-25, init=False, repr=False)
    j: InitVar[float | None] = None  # for spin orbit coupling
    _j: float = field(init=False, repr=False, default=0.0)
    label: InitVar[str | None] = None
    _label: str = field(init=False, repr=False, default="")
    alpha: str = "UPF"

    def __post_init__(self, j, label):
        """Run the `j` and `label` values through their property setters."""
        self.j = j
        self.label = label

    @property
    def j(self):  # noqa: F811 -- InitVar `j` above feeds this property via __post_init__
        """Total angular momentum number."""
        if np.abs(self._j) < 1e-8:
            raise AttributeError("j is not an attribute.")

        jout = self._j
        return jout

    @j.setter
    def j(self, value):
        if not isinstance(value, float):
            self._j = 0.0
        elif np.abs(self.l - value) - 0.5 < 1e-8:
            self._j = value
        else:
            raise ValueError(f"l={self.l}, j={self._j}")

    @property
    def label(self):  # noqa: F811 -- InitVar `label` above feeds this property via __post_init__
        """Orbital label.

        n_shell+orbital_label(lower case)
        """
        return self._label

    @label.setter
    def label(self, value):
        if not isinstance(value, str):
            self._label = "0" + num2label[self.l]
            return
        if len(value) != 2:
            raise ValueError(f"label<{value}> must be 2-digit")
        nshell = value[0]
        orbital = value[1]
        if not nshell.isdigit():
            raise ValueError(f"first digit of label<{nshell}> must be number")
        if orbital.isalpha():
            orbital = orbital.lower()
            if orbital not in ["s", "p", "d", "f"]:
                raise ValueError("second digit must be s, p, d or f")
        else:
            raise ValueError(f"second digit of label<{nshell}> must be alphabet")

        self._label = nshell + orbital

    # def __init__(self, x, y, l, j=0.0) -> None:
    #     super().__init__(x, y, l)
    #     self.j = j


class newProjectors(Projectors):
    """A list of projectors, with more extra functionalities."""

    def to_str_soc(self) -> str:
        """Convert the Projectors into a string following the format for ``pw2wannier90`` and ``Wannier90``."""
        lines = [
            f"{len(self.data[0].x)} {len(self.data)}",
            " ".join([f"{proj.l:3d}" for proj in self.data]),
            " ".join([f"{proj.j:.1f}" for proj in self.data]),
        ]
        content = np.concatenate(
            [np.vstack([self.data[0].x, self.data[0].r]), [p.y for p in self.data]]
        ).transpose()
        lines += [" ".join([f"{v:18.12e}" for v in row]) for row in content]
        return "\n".join(lines)

    def to_file(self, filename: Path):
        """Dump the Projectors to a file following the format for ``pw2wannier90`` and ``Wannier90``."""
        with open(filename, "w", encoding="utf-8") as fd:
            if hasattr(self[0], "j"):
                fd.write(self.to_str_soc())
            else:
                fd.write(self.to_str())

    @classmethod
    def from_pao(cls, filename: Path | str, n: int, l: int, label: str = None):
        """Create a Projectors object from an openMX flavor PAO file."""

        filename = filename if isinstance(filename, Path) else Path(filename)

        # copy from gen_pao.py
        with open(filename, encoding="utf-8") as fin:
            data = []
            in_tag = False
            for line in fin:
                if f"<pseudo.atomic.orbitals.L={l}" in line:
                    in_tag = True
                elif f"pseudo.atomic.orbitals.L={l}>" in line:
                    in_tag = False
                elif in_tag:
                    data.append(line)

        # have read the data from <pseudo.atomic.orbitals.L={l}>
        x = []
        y = []
        for line in data:
            line_pao = np.fromstring(line, sep=" ")
            # line_pao:
            # x  r  n=1 n=2 n=3 ...
            x.append(line_pao[0])
            y.append(line_pao[n + 2])
        x = np.array(x)
        y = np.array(y)
        data_pao = [newProjector(x, y, l, label=label)]
        return cls(data_pao)

    @classmethod
    def from_str_soc(cls, string: str) -> "Projectors":
        """Create a Projectors object from a string that follows the format for ``pw2wannier90`` and ``Wannier90``."""
        lines = [l for l in string.split("\n") if l]
        content = np.array(
            [[float(v) for v in row.split()] for row in lines[3:]]
        ).transpose()
        lvals = [int(l) for l in lines[1].split()]
        jvals = [float(j) for j in lines[2].split()]
        data = [
            newProjector(content[0], y, l, j)
            for l, j, y in zip(lvals, jvals, content[2:], strict=True)
        ]

        return cls(data)

    def add_projector(self, projector: newProjector):
        """Add a projector element to self."""

        self += [projector]

    def add_projector_soc(self, projector: newProjector):
        """Add a projector element to self.

        If the projector has only l, split it to different j with same radials funtction.
        """
        if hasattr(projector, "j"):
            self += [projector]
        elif projector.l == 0:
            self += [
                newProjector(
                    x=projector.x,
                    y=projector.y,
                    l=projector.l,
                    j=projector.l + 0.5,
                    label=projector.label,
                    alpha=projector.alpha,
                )
            ]
        else:
            projector_minus = newProjector(
                x=projector.x,
                y=projector.y,
                l=projector.l,
                j=projector.l - 0.5,
                label=projector.label,
                alpha=projector.alpha,
            )
            projector_plus = newProjector(
                x=projector.x,
                y=projector.y,
                l=projector.l,
                j=projector.l + 0.5,
                label=projector.label,
                alpha=projector.alpha,
            )
            self += [projector_minus, projector_plus]
