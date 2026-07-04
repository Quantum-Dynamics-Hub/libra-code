from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class EschLevineLinearModel(AnalyticalHamiltonianModel):
    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"nstates": 2, "V": [[0.0, 0.0], [0.0, 0.0]], "w": [[0.0, 0.0], [0.0, 0.0]]})
        x = Q[0]
        V = np.asarray(p["V"])
        w = np.asarray(p["w"])
        n = int(p["nstates"])
        self.nstates = n
        H = zeros_from(x, (*tuple(x.shape), n, n))
        dH = zeros_from(x, (*tuple(x.shape), 1, n, n))
        for i in range(n):
            for j in range(n):
                H[..., i, j] = V[i, j] + w[i, j] * x
                dH[..., 0, i, j] = w[i, j]
        return H, dH


@dataclass
class EschLevineJCP2020Model(AnalyticalHamiltonianModel):
    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 2, "nstates": 2, "delta": 0.01})
        x = Q[0]
        n = int(p["nstates"])
        self.nstates = n
        H = zeros_from(x, (*tuple(x.shape), n, n))
        dH = zeros_from(x, (*tuple(x.shape), 1, n, n))
        H[..., 0, 0] = -p["w0"] * x
        dH[..., 0, 0, 0] = -p["w0"]
        for i in range(1, n):
            shift = p["eps"] if i >= p["i_crit"] else 0.0
            H[..., i, i] = p["w1"] * x - i * p["delta"] - shift
            H[..., 0, i] = p["V"]
            H[..., i, 0] = p["V"]
            dH[..., 0, i, i] = p["w1"]
        return H, dH


def esch_levine_jcp2020_params(set_index):
    if set_index == 1:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 9, "delta": 0.01}
    if set_index == 2:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.08, "i_crit": 5, "nstates": 9, "delta": 0.01}
    if set_index == 3:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.005}
    if set_index == 4:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.0025}
    if set_index == 5:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.00125}
    if set_index == 6:
        return {"w0": 0.25, "w1": 0.025, "V": 0.005, "eps": 0.0, "i_crit": 0, "nstates": 17, "delta": 0.000625}
    raise ValueError("set_index must be in 1..6")
