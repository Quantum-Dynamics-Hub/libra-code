from __future__ import annotations

from dataclasses import dataclass

from ._helpers import zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class LVCModel(AnalyticalHamiltonianModel):
    nstates: int = 2

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        omega = params["omega"]
        d1 = params["d1"]
        d2 = params["d2"]
        coup = params["coup"]
        mass = params["mass"]
        ndof = len(omega)
        self.ndof = ndof
        H = zeros_from(Q[0], (*tuple(Q[0].shape), 2, 2))
        dH = zeros_from(Q[0], (*tuple(Q[0].shape), ndof, 2, 2))
        H[..., 0, 0] = params["Delta1"]
        H[..., 1, 1] = params["Delta2"]
        for n in range(ndof):
            qn = Q[n]
            sqrt_m = xp.sqrt(xp.asarray(mass[n]))
            bath = 0.5 * mass[n] * omega[n] * omega[n] * qn * qn
            dbath = mass[n] * omega[n] * omega[n] * qn
            H[..., 0, 0] = H[..., 0, 0] + bath + sqrt_m * qn * d1[n]
            H[..., 1, 1] = H[..., 1, 1] + bath + sqrt_m * qn * d2[n]
            H[..., 0, 1] = H[..., 0, 1] + sqrt_m * coup[n] * qn
            H[..., 1, 0] = H[..., 0, 1]
            dH[..., n, 0, 0] = dbath + sqrt_m * d1[n]
            dH[..., n, 1, 1] = dbath + sqrt_m * d2[n]
            dH[..., n, 0, 1] = sqrt_m * coup[n]
            dH[..., n, 1, 0] = sqrt_m * coup[n]
        return H, dH
