from __future__ import annotations

from dataclasses import dataclass

from .base import AnalyticalHamiltonianModel
from ._helpers import two_state_with_derivatives


@dataclass
class TullyModel1(AnalyticalHamiltonianModel):
    """Tully simple avoided crossing model."""

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
    """Tully dual avoided crossing model."""

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
    """Tully extended coupling with reflection model."""

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
