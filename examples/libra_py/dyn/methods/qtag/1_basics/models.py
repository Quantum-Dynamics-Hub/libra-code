# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Pure-array versions of the tutorial's Holstein and Ferretti models.

Existing dyn.models provide values/first derivatives. Here we add analytic
Hessians so the same example can also exercise LHA potential quadrature.
Electronic axes are last; all arrays stay on the coordinate backend/device.
"""
import numpy as np
from libra_py.dyn.models import Holstein5Model, FerrettiModel
from libra_py.dyn.methods.qtag.core import HamiltonianData


def compute_model(q, params):
    """Return H(n,2,2), dH(n,d,2,2), d2H(n,d,d,2,2)."""
    backend = params.get("backend", "numpy")
    if backend == "torch":
        import torch as xp
        def zero(shape):
            return xp.zeros(shape, dtype=q.dtype, device=q.device)
    else:
        xp = np
        def zero(shape):
            return np.zeros(shape, dtype=q.dtype)
    if params.get("case", "holstein") == "holstein":
        model = Holstein5Model(params=params, backend=backend, nstates=2)
        H, dH = model.diabatic_with_derivatives(q, params)
        hess = zero((q.shape[1], 1, 1, 2, 2))
        for state in range(2):
            hess[:, 0, 0, state, state] = params["k_n"][state]
        for i, j in ((0, 1), (1, 0)):
            alpha = params["alpha"][i][j]
            dx = q[0]-params["x_nm"][i][j]
            hess[:, 0, 0, i, j] = (4*alpha**2*dx**2-2*alpha)*H[:, i, j].real
    else:
        model = FerrettiModel(params=params, backend=backend)
        H, dH = model.diabatic_with_derivatives(q, params)
        hess = zero((q.shape[1], 2, 2, 2, 2))
        for state in range(2):
            hess[:, 0, 0, state, state] = params["Kx"]
            hess[:, 1, 1, state, state] = params["Ky"]
        dx, y = q[0]-params["X3"], q[1]
        alpha, beta, gamma = params["alpha"], params["beta"], params["gamma"]
        expxy = xp.exp(-alpha*dx**2-beta*y**2)
        for i, j in ((0, 1), (1, 0)):
            hess[:, 0, 0, i, j] = (4*alpha**2*dx**2-2*alpha)*H[:, i, j].real
            hess[:, 1, 1, i, j] = gamma*expxy*(-6*beta*y+4*beta**2*y**3)
            hess[:, 0, 1, i, j] = -2*alpha*dx*gamma*expxy*(1-2*beta*y**2)
            hess[:, 1, 0, i, j] = hess[:, 0, 1, i, j]
    return HamiltonianData(H, dH, hess)
