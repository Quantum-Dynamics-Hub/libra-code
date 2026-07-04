from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class Holstein2Model(AnalyticalHamiltonianModel):
    """
    Generic Holstein model with nearest-neighbor constant electronic coupling.

    ``Q[0]=x``. Each diabatic state is a displaced harmonic surface
    ``E_n + 1/2 k_n (x-x_n)^2``. Neighboring states are coupled by the constant
    ``V``. State count and surfaces come from ``E_n``, ``x_n``, and ``k_n``.

    Legacy source: ``libra_py.models.Holstein.Holstein2``.
    """

    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"V": 0.001})
        return _holstein_parabolic(Q[0], self, p, coupling="nearest")


@dataclass
class Holstein3Model(AnalyticalHamiltonianModel):
    """
    Generic Holstein model with distance-indexed constant couplings.

    Diagonal surfaces are displaced harmonic potentials. Off-diagonal coupling
    between states ``i`` and ``j`` is selected from ``V_n[abs(i-j)-1]``.

    Legacy source: ``libra_py.models.Holstein.Holstein3``.
    """

    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {"V_n": [0.001, 0.0, 0.0]})
        return _holstein_parabolic(Q[0], self, p, coupling="distance")


@dataclass
class Holstein4Model(AnalyticalHamiltonianModel):
    """
    Generic Holstein model with a full coupling matrix.

    Diagonal surfaces are displaced harmonic potentials. Off-diagonal elements
    are read directly from ``V[i][j]``. This is useful when couplings do not
    depend only on state-index distance.

    Legacy source: ``libra_py.models.Holstein.Holstein4``.
    """

    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        p = merged_params(params, {})
        return _holstein_parabolic(Q[0], self, p, coupling="matrix")


@dataclass
class Holstein5Model(AnalyticalHamiltonianModel):
    """
    Holstein model with Gaussian off-diagonal couplings.

    Diagonal surfaces are displaced harmonic potentials
    ``E_n + 1/2 k_n (x-x_n)^2``. Off-diagonal terms are
    ``V_ij exp[-alpha_ij (x-x_nm,ij)^2]`` with analytical derivatives.

    Legacy source: ``libra_py.models.Holstein.Holstein5``.
    """

    ndof: int = 1

    def diabatic_with_derivatives(self, Q, params):
        xp = self.xp
        p = merged_params(params, {})
        x = Q[0]
        E, x_n, k_n = p["E_n"], p["x_n"], p["k_n"]
        V, alpha, x_nm = p["V"], p["alpha"], p["x_nm"]
        n = len(E)
        self.nstates = n
        H = zeros_from(x, (*tuple(x.shape), n, n))
        dH = zeros_from(x, (*tuple(x.shape), 1, n, n))
        for i in range(n):
            H[..., i, i] = E[i] + 0.5 * k_n[i] * (x - x_n[i]) ** 2
            dH[..., 0, i, i] = k_n[i] * (x - x_n[i])
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                z = V[i][j] * xp.exp(-alpha[i][j] * (x - x_nm[i][j]) ** 2)
                H[..., i, j] = z
                dH[..., 0, i, j] = -2.0 * alpha[i][j] * (x - x_nm[i][j]) * z
        return H, dH


def _holstein_parabolic(x, model, params, coupling):
    E, x_n, k_n = params["E_n"], params["x_n"], params["k_n"]
    n = len(E)
    model.nstates = n
    H = zeros_from(x, (*tuple(x.shape), n, n))
    dH = zeros_from(x, (*tuple(x.shape), 1, n, n))
    for i in range(n):
        H[..., i, i] = E[i] + 0.5 * k_n[i] * (x - x_n[i]) ** 2
        dH[..., 0, i, i] = k_n[i] * (x - x_n[i])

    if coupling == "nearest":
        for i in range(n - 1):
            H[..., i, i + 1] = params["V"]
            H[..., i + 1, i] = params["V"]
    elif coupling == "distance":
        V_n = params["V_n"]
        for i in range(n):
            for j in range(n):
                if i != j:
                    H[..., i, j] = V_n[abs(i - j) - 1]
    elif coupling == "matrix":
        V = np.asarray(params["V"])
        for i in range(n):
            for j in range(n):
                if i != j:
                    H[..., i, j] = V[i, j]
    return H, dH
