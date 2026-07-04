from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, two_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class LibraModel1(AnalyticalHamiltonianModel):
    """
    Libra internal spin-boson/Marcus two-state model.

    ``Q[0]=x``. The diabatic matrix is
    ``[[k x^2, V], [V, k(x-x0)^2 + D]]`` with analytical derivatives.
    ``x0`` is the displacement between diabatic minima, ``k`` is the force
    constant, ``D`` is the energy bias, and ``V`` is the electronic coupling.

    Legacy source: ``libra_py.models.Libra.model1``. The same form is related
    in the legacy module to the Landry-Subotnik spin-boson parameterization.
    """

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"x0": 1.0, "k": 0.01, "D": 0.0, "V": 0.005})
        x = Q[0]
        h00 = p["k"] * x * x
        h11 = p["k"] * (x - p["x0"]) ** 2 + p["D"]
        h01 = xp.zeros_like(x) + p["V"]
        dh00 = 2.0 * p["k"] * x
        dh11 = 2.0 * p["k"] * (x - p["x0"])
        dh01 = xp.zeros_like(x)
        return two_state_with_derivatives(self, x, h00, h11, h01, dh00, dh11, dh01)
