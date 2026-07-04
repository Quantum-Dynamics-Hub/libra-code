from __future__ import annotations

from dataclasses import dataclass

from .base import AnalyticalHamiltonianModel
from ._helpers import two_state_with_derivatives


@dataclass
class TullyModel1(AnalyticalHamiltonianModel):
    """
    Tully model 1: simple avoided crossing (SAC).

    ``Q[0] = x``. The diabatic Hamiltonian is
    ``H00 = A(1-exp(-Bx))`` for ``x>0`` and ``H00 = -A(1-exp(Bx))`` for
    ``x<=0``, ``H11 = -H00``, and ``H01 = C exp(-D x^2)``. Derivatives are the
    analytical ``dH/dx``. Defaults are the canonical Tully parameters in
    atomic units.

    Reference: J. C. Tully, J. Chem. Phys. 1990, 93, 1061.
    Legacy source: ``libra_py.models.Tully.Tully1_py``.
    """

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        x = Q[0]
        A = params.get("A", 0.010)
        B = params.get("B", 1.600)
        C = params.get("C", 0.005)
        D = params.get("D", 1.000)

        exp_pos = xp.exp(-B * x)
        exp_neg = xp.exp(B * x)
        V11 = xp.where(x > 0, A * (1.0 - exp_pos), -A * (1.0 - exp_neg))
        dV11 = xp.where(x > 0, A * B * exp_pos, A * B * exp_neg)

        exp_c = xp.exp(-D * x * x)
        V12 = C * exp_c
        dV12 = -2.0 * x * C * D * exp_c

        return two_state_with_derivatives(self, x, V11, -V11, V12, dV11, -dV11, dV12)


@dataclass
class TullyModel2(AnalyticalHamiltonianModel):
    """
    Tully model 2: dual avoided crossing (DAC).

    ``Q[0] = x``. The diabatic Hamiltonian is ``H00 = 0``,
    ``H11 = E - A exp(-B x^2)``, and ``H01 = C exp(-D x^2)``. The parameter
    name ``E0`` is accepted as an alias for ``E``. Derivatives are analytical.

    Reference: J. C. Tully, J. Chem. Phys. 1990, 93, 1061.
    Legacy source: ``libra_py.models.Tully.Tully2``.
    """

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        x = Q[0]
        A = params.get("A", 0.100)
        B = params.get("B", 0.280)
        C = params.get("C", 0.015)
        D = params.get("D", 0.060)
        E = params.get("E", params.get("E0", 0.050))

        z_h = xp.exp(-B * x * x)
        z_c = xp.exp(-D * x * x)
        V11 = xp.zeros_like(x)
        V22 = E - A * z_h
        V12 = C * z_c
        dV11 = xp.zeros_like(x)
        dV22 = 2.0 * A * B * x * z_h
        dV12 = -2.0 * C * D * x * z_c

        return two_state_with_derivatives(self, x, V11, V22, V12, dV11, dV22, dV12)


@dataclass
class TullyModel3(AnalyticalHamiltonianModel):
    """
    Tully model 3: extended coupling with reflection (ECWR).

    ``Q[0] = x``. The diagonal elements are constant, ``H00=A`` and
    ``H11=-A``. The off-diagonal coupling is exponential on the left,
    ``B exp(Cx)``, and approaches ``2B`` on the right as
    ``B(2-exp(-Cx))``. The derivative tensor stores ``dH01/dx``.

    Reference: J. C. Tully, J. Chem. Phys. 1990, 93, 1061.
    Legacy source: ``libra_py.models.Tully.Tully3``.
    """

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        x = Q[0]
        A = params.get("A", 0.0006)
        B = params.get("B", 0.1000)
        C = params.get("C", 0.9000)

        V11 = xp.zeros_like(x) + A
        V22 = -V11
        exp_left = xp.exp(C * x)
        exp_right = xp.exp(-C * x)
        V12 = xp.where(x <= 0, B * exp_left, B * (2.0 - exp_right))
        dV12 = xp.where(x <= 0, B * C * exp_left, B * C * exp_right)
        dV11 = xp.zeros_like(x)
        dV22 = xp.zeros_like(x)

        return two_state_with_derivatives(self, x, V11, V22, V12, dV11, dV22, dV12)
