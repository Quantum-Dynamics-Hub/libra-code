from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, two_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class SubotnikDumbbellModel(AnalyticalHamiltonianModel):
    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"A": 0.0006, "B": 0.1, "C": 0.9, "Z": 10.0})
        x = Q[0]
        B, C, Z = p["B"], p["C"], p["Z"]
        left = x < -Z
        middle = (x >= -Z) & (x <= Z)
        e_pp = xp.exp(C * (x + Z))
        e_pm = xp.exp(C * (x - Z))
        e_np = xp.exp(-C * (x + Z))
        e_nm = xp.exp(-C * (x - Z))
        h01 = xp.where(left, B * (e_pm + 2.0 - e_pp), xp.where(middle, B * (e_pm + e_np), B * (e_np + 2.0 - e_nm)))
        dh01 = xp.where(left, B * C * (e_pm - e_pp), xp.where(middle, B * C * (e_pm - e_np), -B * C * (e_np - e_nm)))
        zero = xp.zeros_like(x)
        return two_state_with_derivatives(self, x, zero + p["A"], zero - p["A"], h01, zero, zero, dh01)


@dataclass
class SubotnikDoubleArchModel(AnalyticalHamiltonianModel):
    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {"A": 0.0006, "B": 0.1, "C": 0.9, "Z": 4.0})
        x = Q[0]
        B, C, Z = p["B"], p["C"], p["Z"]
        left = x < -Z
        middle = (x >= -Z) & (x <= Z)
        e_pp = xp.exp(C * (x + Z))
        e_pm = xp.exp(C * (x - Z))
        e_np = xp.exp(-C * (x + Z))
        e_nm = xp.exp(-C * (x - Z))
        h01 = xp.where(left, B * (-e_pm + e_pp), xp.where(middle, B * (-e_pm - e_np + 2.0), B * (e_nm - e_np)))
        dh01 = xp.where(left, B * C * (-e_pm + e_pp), xp.where(middle, B * C * (-e_pm + e_np), B * C * (-e_nm + e_np)))
        zero = xp.zeros_like(x)
        return two_state_with_derivatives(self, x, zero + p["A"], zero - p["A"], h01, zero, zero, dh01)
