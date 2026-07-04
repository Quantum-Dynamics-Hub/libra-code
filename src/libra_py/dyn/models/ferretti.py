from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class FerrettiModel(AnalyticalHamiltonianModel):
    ndof: int = 2
    nstates: int = 2

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(
            params,
            {"X1": 4.0, "X2": 3.0, "X3": 3.0, "Kx": 0.02, "Ky": 0.1, "Delta": 0.01, "alpha": 3.0, "beta": 1.5, "gamma": 0.005},
        )
        X, Y = Q[0], Q[1]
        H = zeros_from(X, (*tuple(X.shape), 2, 2))
        dH = zeros_from(X, (*tuple(X.shape), 2, 2, 2))
        expx = xp.exp(-p["alpha"] * (X - p["X3"]) ** 2)
        expy = xp.exp(-p["beta"] * Y * Y)
        h12 = p["gamma"] * Y * expx * expy
        H[..., 0, 0] = 0.5 * p["Kx"] * (X - p["X1"]) ** 2 + 0.5 * p["Ky"] * Y * Y
        H[..., 1, 1] = 0.5 * p["Kx"] * (X - p["X2"]) ** 2 + 0.5 * p["Ky"] * Y * Y + p["Delta"]
        H[..., 0, 1] = h12
        H[..., 1, 0] = h12
        dh12_x = -2.0 * p["alpha"] * (X - p["X3"]) * h12
        dh12_y = p["gamma"] * expx * expy * (1.0 - 2.0 * p["beta"] * Y * Y)
        dH[..., 0, 0, 0] = p["Kx"] * (X - p["X1"])
        dH[..., 0, 1, 1] = p["Kx"] * (X - p["X2"])
        dH[..., 0, 0, 1] = dh12_x
        dH[..., 0, 1, 0] = dh12_x
        dH[..., 1, 0, 0] = p["Ky"] * Y
        dH[..., 1, 1, 1] = p["Ky"] * Y
        dH[..., 1, 0, 1] = dh12_y
        dH[..., 1, 1, 0] = dh12_y
        return H, dH
