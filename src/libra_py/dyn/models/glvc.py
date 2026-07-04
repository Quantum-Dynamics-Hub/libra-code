from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class GLVCModel(AnalyticalHamiltonianModel):
    """Generalized linear vibronic coupling model."""

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"Ham": [[0.0, 0.0], [0.0, 0.0]], "nstates": 2, "num_osc": 2, "coupling_scaling": [1.0, -1.0], "mass": [1.0, 1.0]})
        w = p["omega"]
        coupl = p["coupl"]
        nstates = int(p["nstates"])
        num_osc = int(p["num_osc"])
        Ham = np.asarray(p["Ham"])
        scale = p["coupling_scaling"]
        mass = p["mass"]
        self.ndof = num_osc
        self.nstates = nstates
        H = zeros_from(Q[0], (*tuple(Q[0].shape), nstates, nstates))
        dH = zeros_from(Q[0], (*tuple(Q[0].shape), num_osc, nstates, nstates))
        for i in range(nstates):
            for j in range(nstates):
                H[..., i, j] = Ham[i, j]
        for n in range(nstates):
            for f in range(num_osc):
                qf = Q[f]
                w2 = mass[f] * w[n][f] ** 2
                H[..., n, n] = H[..., n, n] + 0.5 * w2 * qf * qf + coupl[n][f] * qf * scale[n]
                dH[..., f, n, n] = dH[..., f, n, n] + w2 * qf + coupl[n][f] * scale[n]
        return H, dH


def glvc_debye_bath_params(nstates=2, num_osc=2, omega_c=0.001, reorganization_energy=0.001):
    omega = [
        omega_c * np.tan(0.5 * np.pi * (1.0 - (k / (num_osc + 1))))
        for k in range(1, num_osc + 1)
    ]
    pref = -np.sqrt(2.0 * reorganization_energy / (num_osc + 1.0))
    coupl = [w * pref for w in omega]
    return {
        "nstates": nstates,
        "num_osc": num_osc,
        "omega": [list(omega) for _ in range(nstates)],
        "coupl": [list(coupl) for _ in range(nstates)],
        "mass": [1.0 for _ in range(num_osc)],
    }
