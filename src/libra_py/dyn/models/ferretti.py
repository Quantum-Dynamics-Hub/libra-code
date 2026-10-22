from __future__ import annotations

from dataclasses import dataclass

from ._helpers import merged_params, zeros_from
from .base import AnalyticalHamiltonianModel


@dataclass
class FerrettiModel(AnalyticalHamiltonianModel):
    """
    Ferretti-Granucci-Lami-Persico-Villani two-state conical-intersection model.

    Coordinates are ``Q[0]=X`` and ``Q[1]=Y``. The diagonal elements are
    displaced harmonic wells in ``X`` with a common harmonic ``Y`` term:
    ``H00 = 1/2 Kx (X-X1)^2 + 1/2 Ky Y^2`` and
    ``H11 = 1/2 Kx (X-X2)^2 + 1/2 Ky Y^2 + Delta``. The coupling is
    ``H01 = gamma Y exp[-alpha (X-X3)^2] exp[-beta Y^2]``, odd in ``Y`` and
    localized in ``X``.

    Reference: A. Ferretti, G. Granucci, A. Lami, M. Persico, and G. Villani,
    J. Chem. Phys. 1996, 104, 5517-5527, https://doi.org/10.1063/1.471791.
    Legacy source: ``libra_py.models.Ferretti``.
    """

    ndof: int = 2
    nstates: int = 2

    def diabatic_hessian(self, Q, params):
        """Analytic Hessian with both mixed X/Y entries; finite at Y=0."""
        p = merged_params(params, {'X3': 3., 'Kx': .02, 'Ky': .1,
                                  'alpha': 3., 'beta': 1.5, 'gamma': .005})
        dx, y = Q[0]-p['X3'], Q[1]
        alpha, beta, gamma = p['alpha'], p['beta'], p['gamma']
        e = self.xp.exp(-alpha*dx**2-beta*y**2)
        out = zeros_from(y, (*tuple(y.shape), 2, 2, 2, 2))
        for state in range(2):
            out[..., 0, 0, state, state] = p['Kx']
            out[..., 1, 1, state, state] = p['Ky']
        for i, j in ((0, 1), (1, 0)):
            out[..., 0, 0, i, j] = (4*alpha**2*dx**2-2*alpha)*gamma*y*e
            out[..., 1, 1, i, j] = gamma*e*(-6*beta*y+4*beta**2*y**3)
            out[..., 0, 1, i, j] = -2*alpha*dx*gamma*e*(1-2*beta*y**2)
            out[..., 1, 0, i, j] = out[..., 0, 1, i, j]
        return out

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
