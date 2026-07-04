from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, one_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class HenonHeilesModel(AnalyticalHamiltonianModel):
    ndof: int = 2
    nstates: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"lam": 0.2})
        x, y = Q[0], Q[1]
        lam = p["lam"]
        lam2 = lam * lam
        r2 = x * x + y * y
        value = 0.5 * r2 + lam * (x * y * y - x * x * x / 3.0) + lam2 * r2 * r2 / 16.0
        d_x = x + lam * (y * y - x * x) + 0.25 * lam2 * r2 * x
        d_y = y + 2.0 * lam * x * y + 0.25 * lam2 * r2 * y
        return one_state_with_derivatives(self, value, (d_x, d_y))
