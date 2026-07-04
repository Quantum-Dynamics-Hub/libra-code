from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class SSYModel(AnalyticalHamiltonianModel):
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
