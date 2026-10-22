# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Electronic evolution of coherent nonorthogonal basis coefficients."""
from ..backends.arrays import asarray, adjoint
from ..utils.linalg import _orthogonalizer
from .electronic import hermitian_propagator

def propagate_coefficients(dt, coeff, H, S, rcond=1e-12):
    """Return exp(-i*S^-1*H*dt)*coeff via a Hermitian generalized eigensolve.

    Small overlap eigenvalues are discarded, not inverted. At fixed rank
    the tensor operations preserve autograd (degenerate eigenspaces may have
    ill-defined eigenvector derivatives). No implicit norm rescaling occurs.
    """
    c = asarray(coeff, S, True)
    shape = c.shape
    c = c.reshape(S.shape[0], -1)
    W = _orthogonalizer(S, rcond)
    h = adjoint(W) @ H @ W
    return hermitian_propagator(c, h, dt, basis=W, metric=S).reshape(shape)
