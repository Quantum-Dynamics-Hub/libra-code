# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Normalized Gaussian expansions, using alpha/2 in the exponent.

Kernel arrays use (ndof, npackets), matching dyn.gwp.integrals. Storage-facing
code transposes (npackets, ndof) slices at this boundary. Width alpha is twice
the GWP integral width; phase is per DOF. All calculations are in atomic units.
"""
import numpy as np
from . import integrals as gwp
from ..backends.arrays import xp_of, asarray

def validate_basis(q, p, alp, s):
    xp, a = gwp._arrays((q, p, alp, s))
    q, p, alp, s = a
    if q.ndim != 2 or 0 in q.shape or any(x.shape != q.shape for x in a):
        raise ValueError("q,p,alp,s must share a nonempty (ndof,ntraj) shape")
    if not all(bool(xp.isfinite(x).all()) for x in a):
        raise ValueError("packet parameters must be finite")
    gwp._positive(alp, "alp")
    return a


def validate_states(states, ntraj, nstates=None):
    vals = list(states)
    if not vals or len(vals) != ntraj or any(int(x) != x or x < 0 for x in vals):
        raise ValueError("one nonnegative integer electronic label is required per packet")
    vals = [int(x) for x in vals]
    if nstates is not None and max(vals) >= nstates:
        raise ValueError("electronic label exceeds model state count")
    return vals


def evaluate_packets(points, q, p, alp, s):
    """Return G_j(r_k), shape (ntraj,npoints), and displacements (d,n,m)."""
    points = asarray(points, q)
    if points.ndim == 1:
        points = points[:, None]
    if points.ndim != 2 or points.shape[0] != q.shape[0]:
        raise ValueError("points must have shape (ndof,npoints) or (ndof,)")
    xp = xp_of(q)
    dq = points[:, None, :] - q[:, :, None]
    log_g = (0.25*xp.log(alp[:, :, None]/np.pi)-0.5*alp[:, :, None]*dq**2
             + 1j*(p[:, :, None]*dq+s[:, :, None]))
    return xp.exp(log_g.sum(axis=0)), dq


def wavefunction(points, q, p, alp, s, Coeff):
    """Evaluate wavefunctions at one or many points.

    Coeff is (nstates,npackets); output is (nstates,npoints). For state-labeled
    packets, set entries on other states to zero. All phases are included,
    including during saved reconstruction.
    """
    q, p, alp, s = validate_basis(q, p, alp, s)
    c = asarray(Coeff, q, complex=True)
    if c.ndim != 2 or c.shape[1] != q.shape[1]:
        raise ValueError("Coeff must have shape (nstates,ntraj)")
    return c @ evaluate_packets(points, q, p, alp, s)[0]


def spatial_overlap(q, p, alp, s):
    """Return (ntraj,ntraj) spatial overlaps, without electronic labels."""
    q, p, alp, s = validate_basis(q, p, alp, s)
    return gwp.gwp_overlap_matrix(q, p, s, alp/2, q, p, s, alp/2)


def spatial_kinetic(q, p, alp, s, invM):
    """Return spatial kinetic integrals; invM is a diagonal inverse-mass vector."""
    q, p, alp, s = validate_basis(q, p, alp, s)
    return gwp.gwp_kinetic_matrix(q, p, s, alp/2, q, p, s, alp/2, invM)


def momentum_fields(q, p, alp, s, coeff, node_tol=1e-14):
    """Return log-gradient, its diagonal spatial derivative, and density.

    Evaluate at every center: outputs (ndof,ntraj), (ndof,ntraj), (ntraj,).
    log-gradient = grad(psi)/psi, with imaginary part the Bohmian momentum.
    Use the expansion packet's width in its normalization (correcting C++).
    Near nodes (|psi| <= node_tol*max|psi|), ratios are set to zero; fitting
    density weights suppresses these undefined samples. Zero-density surfaces
    should be synchronized/frozen by the workflow rather than propagated.
    """
    q, p, alp, s = validate_basis(q, p, alp, s)
    xp = xp_of(q)
    c = asarray(coeff, q, True).reshape(-1)
    if c.shape[0] != q.shape[1]:
        raise ValueError("one coefficient per packet is required")
    g, dq = evaluate_packets(q, q, p, alp, s)
    weights = c[:, None]*g
    psi = weights.sum(axis=0)
    factors = -alp[:, :, None]*dq+1j*p[:, :, None]
    dpsi = (factors*weights[None]).sum(axis=1)
    d2psi = ((factors**2-alp[:, :, None])*weights[None]).sum(axis=1)
    valid = abs(psi) > node_tol*abs(psi).max()
    safe = xp.where(valid, psi, xp.ones_like(psi))
    ratio = xp.where(valid[None], dpsi/safe[None], xp.zeros_like(dpsi))
    slope = xp.where(valid[None], d2psi/safe[None]-ratio**2, xp.zeros_like(d2psi))
    return ratio, slope, abs(psi)**2
