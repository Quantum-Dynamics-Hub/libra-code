from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class SSYModel(AnalyticalHamiltonianModel):
    """
    Shenvi-Subotnik-Yang two-state, two-dimensional model.

    Coordinates are ``Q[0]=x`` and ``Q[1]=y``. The first diabatic state has
    constant energy ``H00=-E0``. The second has a tilted Gaussian well in the
    rotated coordinates ``x+y`` and ``x-y``. The off-diagonal coupling is
    another anisotropic Gaussian:

    ``H11 = -A exp[-B(0.75(x+y)^2 + 0.25(x-y)^2)]``
    ``H01 = C exp[-D(0.25(x+y)^2 + 0.75(x-y)^2)]``.

    Analytical derivatives with respect to both coordinates are returned.

    Reference: N. Shenvi, J. E. Subotnik, and W. Yang,
    J. Chem. Phys. 2011, 135, 024101. Legacy source:
    ``libra_py.models.SSY.SSY``.
    """

    ndof: int = 2
    nstates: int = 2

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"E0": 0.05, "A": 0.15, "B": 0.14, "C": 0.015, "D": 0.06})
        x, y = Q[0], Q[1]
        H = zeros_from(x, (*tuple(x.shape), 2, 2))
        dH = zeros_from(x, (*tuple(x.shape), 2, 2, 2))
        H[..., 0, 0] = -p["E0"]
        z = -p["A"] * xp.exp(-p["B"] * (0.75 * (x + y) ** 2 + 0.25 * (x - y) ** 2))
        H[..., 1, 1] = z
        dH[..., 0, 1, 1] = -p["B"] * (1.5 * (x + y) + 0.5 * (x - y)) * z
        dH[..., 1, 1, 1] = -p["B"] * (1.5 * (x + y) - 0.5 * (x - y)) * z

        z = p["C"] * xp.exp(-p["D"] * (0.25 * (x + y) ** 2 + 0.75 * (x - y) ** 2))
        dzdx = -p["D"] * (0.5 * (x + y) + 1.5 * (x - y)) * z
        dzdy = -p["D"] * (0.5 * (x + y) - 1.5 * (x - y)) * z
        H[..., 0, 1] = z
        H[..., 1, 0] = z
        dH[..., 0, 0, 1] = dzdx
        dH[..., 0, 1, 0] = dzdx
        dH[..., 1, 0, 1] = dzdy
        dH[..., 1, 1, 0] = dzdy
        return H, dH
