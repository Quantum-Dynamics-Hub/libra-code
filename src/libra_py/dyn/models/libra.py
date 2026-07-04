from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, two_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class LibraModel1(AnalyticalHamiltonianModel):
    """Spin-boson/Marcus two-state model from ``libra_py.models.Libra.model1``."""

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
