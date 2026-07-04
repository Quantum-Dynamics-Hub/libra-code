from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, two_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class ZhuDualRZDModel(AnalyticalHamiltonianModel):
    """
    Zhu dual Rosen-Zener-Demkov model.

    ``Q[0] = x``. The two-state diabatic Hamiltonian uses ``H00=0``,
    ``H11=A``, and a double Gaussian coupling
    ``H01 = B[exp(-C(x-Z)^2) + exp(-C(x+Z)^2)]``. The derivative tensor
    stores the analytical derivative of the coupling.

    Reference: C. Zhu, Sci. Rep. 2016, 6, 24198. Also see
    J. Xu and L. Wang, J. Chem. Phys. 2019, 150, 164101.
    Legacy source: ``libra_py.models.Zhu.dual_RZD``.
    """

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"A": 0.025, "B": 0.025, "C": 0.7, "Z": 3.0})
        x = Q[0]
        e_minus = xp.exp(-p["C"] * (x - p["Z"]) ** 2)
        e_plus = xp.exp(-p["C"] * (x + p["Z"]) ** 2)
        h01 = p["B"] * (e_minus + e_plus)
        dh01 = -2.0 * p["B"] * p["C"] * ((x - p["Z"]) * e_minus + (x + p["Z"]) * e_plus)
        zero = xp.zeros_like(x)
        return two_state_with_derivatives(self, x, zero, zero + p["A"], h01, zero, zero, dh01)


@dataclass
class ZhuDualLZSModel(AnalyticalHamiltonianModel):
    """
    Zhu dual Landau-Zener-Stuckelberg model.

    ``Q[0] = x``. The Hamiltonian is ``H00=0``,
    ``H11 = E0 - A exp(-B x^2)``, and ``H01 = C exp(-D x^2)`` with
    analytical derivatives. It is useful for decoherence and repeated-crossing
    tests.

    Reference: C. Zhu, Sci. Rep. 2016, 6, 24198. Also see
    J. Xu and L. Wang, J. Chem. Phys. 2019, 150, 164101.
    Legacy source: ``libra_py.models.Zhu.dual_LZS``.
    """

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"A": 0.1, "B": 0.28, "C": 0.01, "D": 0.06, "E0": 0.03})
        x = Q[0]
        e_b = xp.exp(-p["B"] * x * x)
        e_d = xp.exp(-p["D"] * x * x)
        h11 = p["E0"] - p["A"] * e_b
        h01 = p["C"] * e_d
        dh11 = 2.0 * p["A"] * p["B"] * e_b * x
        dh01 = -2.0 * p["C"] * p["D"] * e_d * x
        zero = xp.zeros_like(x)
        return two_state_with_derivatives(self, x, zero, h11, h01, zero, dh11, dh01)
