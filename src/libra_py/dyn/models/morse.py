from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class MorseModel(AnalyticalHamiltonianModel):
    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(
            params,
            {
                "D": [0.0, 0.0, 0.0],
                "alpha": [0.0, 0.0, 0.0],
                "x_n": [0.0, 0.0, 0.0],
                "E": [0.0, -0.001, -0.002],
                "V": [[0.001, 0.001, 0.001]] * 3,
                "beta": [[0.0, 0.0, 0.0]] * 3,
                "x_nm": [[0.0, 0.0, 0.0]] * 3,
            },
        )
        x = Q[0]
        E, D, alpha, x_n = p["E"], p["D"], p["alpha"], p["x_n"]
        V, beta, x_nm = np.asarray(p["V"]), p["beta"], p["x_nm"]
        n = len(E)
        self.nstates = n
        H = zeros_from(x, (*tuple(x.shape), n, n))
        dH = zeros_from(x, (*tuple(x.shape), 1, n, n))
        for i in range(n):
            bo = xp.exp(-alpha[i] * (x - x_n[i]))
            H[..., i, i] = E[i] + D[i] * (1.0 - bo) ** 2
            dH[..., 0, i, i] = 2.0 * alpha[i] * D[i] * bo * (1.0 - bo)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                z = V[i, j] * xp.exp(-beta[i][j] * (x - x_nm[i][j]) ** 2)
                H[..., i, j] = z
                dH[..., 0, i, j] = -2.0 * beta[i][j] * (x - x_nm[i][j]) * z
        return H, dH


def coronado_xing_miller_params(model_index):
    params = {}
    if model_index == 1:
        params.update(
            {
                "D": [0.003, 0.004, 0.003],
                "alpha": [0.65, 0.6, 0.65],
                "x_n": [5.00, 4.00, 6.00],
                "E": [0.00, 0.01, 0.006],
                "V": [[0.000, 0.002, 0.000], [0.002, 0.000, 0.002], [0.000, 0.002, 0.000]],
                "x_nm": [[0.00, 3.40, 0.00], [3.40, 0.00, 4.80], [0.00, 4.80, 0.00]],
                "beta": [[0.00, 16.00, 0.00], [16.00, 0.00, 16.00], [0.00, 16.00, 0.00]],
            }
        )
    else:
        raise ValueError("Only Coronado-Xing-Miller model_index=1 is currently provided")
    return params
