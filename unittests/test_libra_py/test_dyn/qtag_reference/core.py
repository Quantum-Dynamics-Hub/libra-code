# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""QTAG Gaussian algebra, potential quadrature and coefficient propagation.

All coordinates use atomic units. Packet arrays q,p,alp,s have shape
(ndof,ntraj). QTAG uses G=(alp/pi)**(1/4)*exp(-alp*(r-q)**2/2+i*p*(r-q)+i*s)
per DOF, so its width is TWICE the width accepted by dyn.gwp.integrals.
Phases s are per DOF. Passing tensors preserves dtype, device and autograd.

The C++ nHamiltonian container is replaced by HamiltonianData or a model
callback returning array data. No C++ extension or SymPy is used.
"""
from dataclasses import dataclass
import numpy as np
from libra_py.dyn.gwp import integrals as gwp
from ._backend import xp_of, asarray, zeros, stack, adjoint, indices

__all__ = ["HamiltonianData", "evaluate_model", "qtag_psi", "qtag_overlap_elementary",
           "qtag_kinetic_elementary", "qtag_overlap", "BAT", "LHA", "BATe", "LHAe",
           "qtag_potential", "qtag_hamiltonian_and_overlap", "qtag_momentum",
           "momentum_fields", "propagate_electronic_qtag", "metric_solve"]


@dataclass
class HamiltonianData:
    """Batched diabatic model values and spatial derivatives.

    H: (ntraj,nstates,nstates); dH: (ntraj,ndof,nstates,nstates);
    d2H: optional (ntraj,ndof,ndof,nstates,nstates). Hessians include BOTH
    mixed-derivative axes. Electronic matrices must be Hermitian. Derivatives
    are derivatives of the potential, not forces. LHA requires d2H.
    """
    H: object
    dH: object
    d2H: object = None


def _basis(q, p, alp, s):
    xp, a = gwp._arrays((q, p, alp, s))
    q, p, alp, s = a
    if q.ndim != 2 or 0 in q.shape or any(x.shape != q.shape for x in a):
        raise ValueError("q,p,alp,s must share a nonempty (ndof,ntraj) shape")
    if not all(bool(xp.isfinite(x).all()) for x in a):
        raise ValueError("packet parameters must be finite")
    gwp._positive(alp, "alp")
    return a


def _states(states, ntraj, nstates=None):
    vals = list(states)
    if len(vals) != ntraj or any(int(x) != x or x < 0 for x in vals):
        raise ValueError("one nonnegative integer electronic label is required per packet")
    vals = [int(x) for x in vals]
    if nstates is not None and max(vals) >= nstates:
        raise ValueError("electronic label exceeds model state count")
    return vals


def evaluate_model(model, q, params=None):
    """Evaluate a model at all packet centers and validate its array contract.

    A callback ``model(q, params)`` returns HamiltonianData, (H,dH[,d2H]),
    or a dict with H_dia/dH_dia[/d2H_dia]. Objects from dyn.models exposing
    evaluate(q, params) work directly for BAT/BATe. No finite-difference
    Hessian is silently substituted when an LHA calculation requests d2H.
    """
    result = model.evaluate(q, params) if hasattr(model, "evaluate") else model(q, params or {})
    if isinstance(result, dict):
        result = HamiltonianData(result["H_dia"], result["dH_dia"], result.get("d2H_dia"))
    elif not isinstance(result, HamiltonianData):
        result = HamiltonianData(*result)
    h = asarray(result.H, q, complex=True)
    dh = asarray(result.dH, q, complex=True)
    d2h = None if result.d2H is None else asarray(result.d2H, q, complex=True)
    n, d = q.shape[1], q.shape[0]
    if h.ndim != 3 or h.shape[0] != n or h.shape[1] != h.shape[2]:
        raise ValueError("H must have shape (ntraj,nstates,nstates)")
    if dh.shape != (n, d, *h.shape[1:]):
        raise ValueError("dH must have shape (ntraj,ndof,nstates,nstates)")
    if d2h is not None and d2h.shape != (n, d, d, *h.shape[1:]):
        raise ValueError("d2H must have shape (ntraj,ndof,ndof,nstates,nstates)")
    xp = xp_of(q)
    for x in (h, dh, d2h):
        if x is not None and (not bool(xp.isfinite(x).all()) or
                              not bool(xp.allclose(x, adjoint(x), rtol=1e-6, atol=1e-8))):
            raise ValueError("model matrices must be finite and Hermitian")
    return HamiltonianData(h, dh, d2h)


def _gaussians(points, q, p, alp, s):
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


def qtag_psi(points, q, p, alp, s, Coeff):
    """Evaluate wavefunctions at one or many points.

    Coeff is (nstates,ntraj); output is (nstates,npoints), matching C++ for
    a single point. For active-state coefficients use coefficient_matrix in
    compute. All phases are included, including during saved reconstruction.
    """
    q, p, alp, s = _basis(q, p, alp, s)
    c = asarray(Coeff, q, complex=True)
    if c.ndim != 2 or c.shape[1] != q.shape[1]:
        raise ValueError("Coeff must have shape (nstates,ntraj)")
    return c @ _gaussians(points, q, p, alp, s)[0]


def qtag_overlap_elementary(q, p, alp, s):
    """Return (ntraj,ntraj) spatial overlaps, without electronic labels."""
    q, p, alp, s = _basis(q, p, alp, s)
    return gwp.gwp_overlap_matrix(q, p, s, alp/2, q, p, s, alp/2)


def qtag_kinetic_elementary(q, p, alp, s, invM):
    """Return spatial kinetic integrals; invM is a diagonal inverse-mass vector."""
    q, p, alp, s = _basis(q, p, alp, s)
    return gwp.gwp_kinetic_matrix(q, p, s, alp/2, q, p, s, alp/2, invM)


def qtag_overlap(active_states, ovlp, nstates):
    """C++-compatible sparse embedding, shape (ntraj*nstates,)*2.

    Only row/column i*nstates+active_states[i] is populated, and only for
    equal electronic labels. Dynamics uses the compact ntraj-square overlap
    from qtag_hamiltonian_and_overlap instead of this mostly zero embedding.
    """
    labels = _states(active_states, ovlp.shape[0], nstates)
    if ovlp.shape != (len(labels), len(labels)):
        raise ValueError("ovlp must be square")
    idx = indices([i*nstates+n for i, n in enumerate(labels)], ovlp)
    st = indices(labels, ovlp)
    out = zeros((len(labels)*nstates,)*2, ovlp, complex=True)
    out[idx[:, None], idx[None, :]] = ovlp*(st[:, None] == st[None, :])
    return out


def _potential_factor(q1, p1, a1, q2, p2, a2, h1, h2, dh1, dh2,
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


def BAT(Ham1, Ham2, dHam1, dHam2, q1, p1, s1, alp1, n1, q2, p2, s2, alp2, n2):
    """Single-pair BAT potential/overlap ratio with the C++ argument order.

    Ham matrices are (nstates,nstates), dHam arrays are (ndof,nstates,nstates).
    Packet inputs are vectors. s1/s2 cancel in the ratio; multiply the result
    by the spatial overlap to obtain the actual potential matrix element.
    """
    return _single_taylor(Ham1, Ham2, dHam1, dHam2, None, None,
                          q1, p1, alp1, n1, q2, p2, alp2, n2)


def LHA(Ham1, Ham2, dHam1, dHam2, d2Ham1, d2Ham2,
        q1, p1, s1, alp1, n1, q2, p2, s2, alp2, n2):
    """Single-pair LHA ratio, including correct mixed-DOF moments.

    Hessians are (ndof,ndof,nstates,nstates), or flattened C++-style
    (ndof**2,nstates,nstates). Other arguments follow BAT. This expression
    is exact for any quadratic potential, including cross-coordinate terms.
    """
    return _single_taylor(Ham1, Ham2, dHam1, dHam2, d2Ham1, d2Ham2,
                          q1, p1, alp1, n1, q2, p2, alp2, n2)


def _single_taylor(h1, h2, dh1, dh2, hh1, hh2, q1, p1, a1, n1, q2, p2, a2, n2):
    _, values = gwp._arrays((q1, p1, a1, q2, p2, a2))
    q1, p1, a1, q2, p2, a2 = [x.reshape(-1, 1) for x in values]
    d = q1.shape[0]
    h1, h2, dh1, dh2 = [asarray(x, q1, complex=True) for x in (h1, h2, dh1, dh2)]
    hh1 = None if hh1 is None else asarray(hh1, q1, True).reshape(d, d, *h1.shape)[..., n1, n2][None]
    hh2 = None if hh2 is None else asarray(hh2, q1, True).reshape(d, d, *h2.shape)[..., n1, n2][None]
    return _potential_factor(q1, p1, a1, q2, p2, a2,
                             h1[n1, n2].reshape(1), h2[n1, n2].reshape(1),
                             dh1[:, n1, n2][None], dh2[:, n1, n2][None], hh1, hh2)[0, 0]


def _exact_coupling(q1, p1, a1, q2, p2, a2, ABC):
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


def qtag_potential(q1, p1, s1, alp1, n1, traj_on_surf_n1,
                   q2, p2, s2, alp2, n2, traj_on_surf_n2, ham, method=0, ABC=(0., 0., 0.)):
    """Return potential/overlap ratios for two state-specific packet blocks.

    The order follows C++, with HamiltonianData replacing nHamiltonian. Its
    batch indices are given by traj_on_surf_n1/2. Methods: 0 BAT, 1 LHA,
    2 LHA plus exact Gaussian off-diagonal coupling, 3 BAT plus exact coupling.
    Exact coupling means SUM_d A*exp(-B*(r_d-C)**2), not a product in nD.
    All results have shape (ntraj1,ntraj2); multiply by spatial overlaps.
    """
    q1, p1, alp1, s1 = _basis(q1, p1, alp1, s1)
    q2, p2, alp2, s2 = _basis(q2, p2, alp2, s2)
    if method not in (0, 1, 2, 3):
        raise ValueError("potential method must be 0, 1, 2 or 3")
    if method in (2, 3) and n1 != n2:
        return _exact_coupling(q1, p1, alp1, q2, p2, alp2, ABC)
    i, j = indices(traj_on_surf_n1, q1), indices(traj_on_surf_n2, q1)
    if len(i) != q1.shape[1] or len(j) != q2.shape[1]:
        raise ValueError("one model batch index per packet is required")
    hh1 = hh2 = None
    if method in (1, 2):
        if ham.d2H is None:
            raise ValueError("LHA requires model Hessians d2H")
        hh1, hh2 = ham.d2H[i, :, :, n1, n2], ham.d2H[j, :, :, n1, n2]
    return _potential_factor(q1, p1, alp1, q2, p2, alp2,
                             ham.H[i, n1, n2], ham.H[j, n1, n2],
                             ham.dH[i, :, n1, n2], ham.dH[j, :, n1, n2], hh1, hh2)


def LHAe(i, j, q1, p1, s1, alp1, n1, q2, p2, s2, alp2, n2, AA, BB, CC, ham):
    """Single-pair LHA/exact-coupling ratio; ham is HamiltonianData."""
    args = [x.reshape(-1, 1) for x in (q1, p1, s1, alp1, q2, p2, s2, alp2)]
    return qtag_potential(*args[:4], n1, [i], *args[4:], n2, [j], ham, 2, (AA, BB, CC))[0, 0]


def BATe(i, j, q1, p1, s1, alp1, n1, q2, p2, s2, alp2, n2, AA, BB, CC, ham):
    """Single-pair BAT/exact-coupling ratio; ham is HamiltonianData."""
    args = [x.reshape(-1, 1) for x in (q1, p1, s1, alp1, q2, p2, s2, alp2)]
    return qtag_potential(*args[:4], n1, [i], *args[4:], n2, [j], ham, 3, (AA, BB, CC))[0, 0]


def qtag_hamiltonian_and_overlap(q, p, alp, s, active_states, invM,
                                model, model_params=None, dyn_params=None):
    """Build compact (S,H) matrices, returned rather than mutated.

    model follows evaluate_model. dyn_params selects qtag_pot_approx_method;
    exact-coupling methods read ex_cpl_A/B/C from model_params. State labels
    need not be contiguous or sorted. Kinetic and overlap blocks couple only
    equal electronic states; potentials couple all states present in the basis.
    """
    q, p, alp, s = _basis(q, p, alp, s)
    data = evaluate_model(model, q, model_params)
    states = _states(active_states, q.shape[1], data.H.shape[1])
    method = (dyn_params or {}).get("qtag_pot_approx_method", 0)
    pars = model_params or {}
    if method in (2, 3) and not all(k in pars for k in ("ex_cpl_A", "ex_cpl_B", "ex_cpl_C")):
        raise ValueError("exact-coupling methods require ex_cpl_A/B/C")
    ABC = tuple(pars.get(k, 0.) for k in ("ex_cpl_A", "ex_cpl_B", "ex_cpl_C"))
    spatial = qtag_overlap_elementary(q, p, alp, s)
    labels = indices(states, q)
    mask = labels[:, None] == labels[None, :]
    S = spatial*mask
    H = qtag_kinetic_elementary(q, p, alp, s, invM)*mask
    for n1 in sorted(set(states)):
        i = indices([k for k, n in enumerate(states) if n == n1], q)
        for n2 in sorted(set(states)):
            j = indices([k for k, n in enumerate(states) if n == n2], q)
            v = qtag_potential(q[:, i], p[:, i], s[:, i], alp[:, i], n1, i,
                               q[:, j], p[:, j], s[:, j], alp[:, j], n2, j, data, method, ABC)
            H[i[:, None], j[None, :]] = H[i[:, None], j[None, :]]+v*spatial[i[:, None], j[None, :]]
    return S, H


def momentum_fields(q, p, alp, s, coeff, node_tol=1e-14):
    """Return log-gradient, its diagonal spatial derivative, and density.

    Evaluate at every center: outputs (ndof,ntraj), (ndof,ntraj), (ntraj,).
    log-gradient = grad(psi)/psi, with imaginary part the Bohmian momentum.
    Use the expansion packet's width in its normalization (correcting C++).
    Near nodes (|psi| <= node_tol*max|psi|), ratios are set to zero; fitting
    density weights suppresses these undefined samples. Zero-density surfaces
    should be synchronized/frozen by the workflow rather than propagated.
    """
    q, p, alp, s = _basis(q, p, alp, s)
    xp = xp_of(q)
    c = asarray(coeff, q, True).reshape(-1)
    if c.shape[0] != q.shape[1]:
        raise ValueError("one coefficient per packet is required")
    g, dq = _gaussians(q, q, p, alp, s)
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


def qtag_momentum(q, p, alp, s, Coeff):
    """Return grad(psi)/psi, shape (ndof,ntraj); momentum is its imaginary part."""
    return momentum_fields(q, p, alp, s, Coeff)[0]


def _orthogonalizer(S, rcond):
    xp = xp_of(S)
    if not 0 <= rcond < 1:
        raise ValueError("overlap_rcond must lie in [0,1)")
    vals, vecs = xp.linalg.eigh((S+adjoint(S))/2)
    keep = vals > rcond*vals.max()
    if not bool(keep.any()):
        raise ValueError("overlap has no positive retained eigenvalues")
    return vecs[:, keep]/xp.sqrt(vals[keep])[None, :]


def metric_solve(S, rhs, rcond=1e-12):
    """Solve S*x=rhs in its positive spectral subspace (relative cutoff)."""
    W = _orthogonalizer(S, rcond)
    return W @ (adjoint(W) @ rhs)


def propagate_electronic_qtag(dt, coeff, H, S, rcond=1e-12):
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
