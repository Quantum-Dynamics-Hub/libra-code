from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, one_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class HenonHeilesModel(AnalyticalHamiltonianModel):
    """
    Two-dimensional Henon-Heiles potential.

    ``Q[0]=x`` and ``Q[1]=y``. The single-state potential is

    ``V = 1/2(x^2+y^2) + lam(x y^2 - x^3/3) + lam^2 (x^2+y^2)^2/16``.

    The returned derivative tensor contains the analytical gradients with
    respect to ``x`` and ``y``. The diabatic overlap is the scalar identity and
    derivative couplings are zero.

    Reference: E. Sim and N. Makri, J. Chem. Phys. 1995, 102, 5616-5625.
    Legacy source: ``libra_py.models.Henon_Heiles.Henon_Heiles``.
    """

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
