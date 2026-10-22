# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Electronic evolution of coherent nonorthogonal basis coefficients."""
from ..backends.arrays import xp_of, asarray, adjoint
from ..utils.linalg import _orthogonalizer

def propagate_coefficients(dt, coeff, H, S, rcond=1e-12):
    """Return exp(-i*S^-1*H*dt)*coeff via a Hermitian generalized eigensolve.

    Small overlap eigenvalues are discarded, not inverted. At fixed rank
    the tensor operations preserve autograd (degenerate eigenspaces may have
    ill-defined eigenvector derivatives). No implicit norm rescaling occurs.
    """
    xp = xp_of(S)
    c = asarray(coeff, S, True)
    shape = c.shape
    c = c.reshape(S.shape[0], -1)
    W = _orthogonalizer(S, rcond)
    h = adjoint(W) @ H @ W
    vals, U = xp.linalg.eigh((h+adjoint(h))/2)
    Z = W @ U
    return (Z @ (xp.exp(-1j*dt*vals)[:, None]*(adjoint(Z) @ S @ c))).reshape(shape)
