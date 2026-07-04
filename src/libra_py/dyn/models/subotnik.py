from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, two_state_with_derivatives
from .base import AnalyticalHamiltonianModel


@dataclass
class SubotnikDumbbellModel(AnalyticalHamiltonianModel):
    """
    Subotnik-Shenvi dumbbell geometry.

    ``Q[0] = x``. This two-state Hamiltonian is a symmetrized version of
    Tully's extended-coupling-with-reflection model. The diagonal elements are
    ``H00=A`` and ``H11=-A``. The coupling is piecewise exponential with two
    centers at ``+-Z`` and scale ``C``; its value and derivative are continuous
    across the central region.

    Parameters ``A``, ``B``, ``C``, and ``Z`` follow the legacy defaults.
    Reference: J. E. Subotnik and N. Shenvi, J. Chem. Phys. 2011,
    134, 024105. Also see J. Xu and L. Wang, J. Chem. Phys. 2019,
    150, 164101.
    """

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
    """
    Subotnik-Shenvi double-arch geometry.

    ``Q[0] = x``. The diagonal terms are ``H00=A`` and ``H11=-A``. The
    off-diagonal element is a piecewise combination of exponentials arranged to
    produce two coupling arches between ``-Z`` and ``Z``. The derivative tensor
    stores the analytical derivative of the coupling with respect to ``x``.

    Reference: J. E. Subotnik and N. Shenvi, J. Chem. Phys. 2011,
    134, 024105. Also see J. Xu and L. Wang, J. Chem. Phys. 2019,
    150, 164101.
    """

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
