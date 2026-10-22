# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Gaussian averages of Taylor potentials and exact Gaussian couplings.

These kernels return potential/overlap ratios, with rows indexing bra packets
and columns ket packets. They use the alpha/2 Gaussian convention and include
both mixed Hessian axes. No storage, model callback, or propagation is owned here.
"""
from ..backends.arrays import xp_of, asarray

def taylor_ratio(q1, p1, a1, q2, p2, a2, h1, h2, dh1, dh2,
                      d2h1=None, d2h2=None):
    """Vectorized bra/ket averaged Taylor expansion divided by overlap."""
    denom = a1[:, :, None]+a2[:, None, :]
    mu = ((a1*q1)[:, :, None]+(a2*q2)[:, None, :]
          + 1j*(p2[:, None, :]-p1[:, :, None]))/denom
    x = mu-q1[:, :, None]
    y = mu-q2[:, None, :]
    v = (h1[:, None]+h2[None, :])/2
    v = v+0.5*(dh1.T[:, :, None]*x+dh2.T[:, None, :]*y).sum(axis=0)
    if d2h1 is not None:
        # E[(r_d-q_d)(r_e-q_e)] = delta_de/(a1_d+a2_d) + x_d*x_e.
        # Sum both mixed axes, so the Taylor factor is 1/2, averaged to 1/4.
        for d in range(q1.shape[0]):
            for e in range(q1.shape[0]):
                variance = 1/denom[d] if d == e else 0
                v = v+0.25*(d2h1[:, d, e][:, None]*(x[d]*x[e]+variance)
                            + d2h2[:, d, e][None, :]*(y[d]*y[e]+variance))
    return v


def gaussian_coupling_ratio(q1, p1, a1, q2, p2, a2, ABC):
    """Ratio for sum_d A*exp[-B*(r_d-C)**2], as in C++ BATe/LHAe."""
    xp = xp_of(q1)
    A, B, C = [asarray(x, q1) for x in ABC]
    if bool((B < 0).any()):
        raise ValueError("Gaussian coupling decay B must be nonnegative")
    denom = a1[:, :, None]+a2[:, None, :]
    mu = ((a1*q1)[:, :, None]+(a2*q2)[:, None, :]
          + 1j*(p2[:, None, :]-p1[:, :, None]))/denom
    scale = 1+2*B/denom
    return (A/xp.sqrt(scale)*xp.exp(-B*(mu-C)**2/scale)).sum(axis=0)
