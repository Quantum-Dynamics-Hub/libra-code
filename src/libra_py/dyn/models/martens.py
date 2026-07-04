from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, one_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class MartensModel1(AnalyticalHamiltonianModel):
    ndof: int = 2
    nstates: int = 1

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"Va": 0.00625, "Vb": 0.0106})
        x, y = Q[0], Q[1]
        sech = 1.0 / xp.cosh(2.0 * x)
        sech2 = sech * sech
        value = p["Va"] * sech2 + 0.5 * p["Vb"] * y * y
        d_x = -4.0 * p["Va"] * xp.tanh(2.0 * x) * sech2
        d_y = p["Vb"] * y
        return one_state_with_derivatives(self, value, (d_x, d_y))


@dataclass
class MartensModel2(AnalyticalHamiltonianModel):
    ndof: int = 2
    nstates: int = 1

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"Va": 0.00625, "Vb": 0.0106, "Vc": 0.4})
        x, y = Q[0], Q[1]
        sech = 1.0 / xp.cosh(2.0 * x)
        sech2 = sech * sech
        shift = y + p["Vc"] * (x * x - 1.0)
        value = p["Va"] * sech2 + 0.5 * p["Vb"] * shift * shift
        d_x = -4.0 * p["Va"] * xp.tanh(2.0 * x) * sech2 + 2.0 * p["Vb"] * p["Vc"] * x * shift
        d_y = p["Vb"] * shift
        return one_state_with_derivatives(self, value, (d_x, d_y))
