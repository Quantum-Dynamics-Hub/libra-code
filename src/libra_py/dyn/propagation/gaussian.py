# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Quantum-momentum fitting and adaptive Gaussian basis evolution.

These functional kernels accept (ndof, npackets) arrays. The QTAG algorithm
reads/writes TensorStorage through the Gaussian view helpers; kernels do not
own trajectory state. Both direct and tutorial-compatible iterative fits are
available. The legacy implementation is preserved separately as a test oracle.
"""
import numpy as np
from ..gwp import integrals as gwp
from ..gwp.packets import validate_basis, momentum_fields
from ..utils.linalg import metric_solve
from ..backends.arrays import xp_of, asarray, stack, copy, indices, scalar

def _legacy_linear_fit(normal, rhs, beta):
    """Reproduce solve_linsys(normal, rhs, x, beta, 200000) in arrays.

    The C++ solver applies Gauss-Seidel to normal.T @ normal, starting
    from zero. It stops when the combined change in the real and imaginary
    fit coefficients, divided by sqrt(2), is at most beta. This tolerance
    changes the moving-basis dynamics; an exact solve is not equivalent.
    """
    if not np.isfinite(beta) or beta <= 0:
        raise ValueError("linfit_beta must be finite and positive")
    xp = xp_of(normal)
    matrix = normal.T @ normal
    target = (normal.T+0j) @ rhs
    if bool((matrix.diagonal() == 0).any()):
        raise ValueError("singular legacy momentum fit; use momentum_fit='direct'")
    fit = xp.zeros_like(rhs)
    for _ in range(200000):
        first = (target[0]-matrix[0, 1]*fit[1])/matrix[0, 0]
        second = (target[1]-matrix[1, 0]*first)/matrix[1, 1]
        updated = stack([first, second])
        error = ((abs(updated-fit)**2).sum()/2)**.5
        fit = updated
        if scalar(error) <= beta:
            return fit
    raise RuntimeError("legacy momentum fit did not converge in 200000 iterations")


def quantum_momentum(q, p, a, s, coeff_on_surf, dyn_params=None):
    """Return momentum, real log-gradient, and their slopes (ndof,ntraj).

    mom_calc_type=0 returns raw ratios (zero slopes, matching legacy updates);
    1 uses a density-weighted linear least-squares fit per coordinate;
    2 returns raw ratios with their analytic spatial derivatives.
    d_weight=0 disables density weighting. momentum_fit='direct' uses a
    spectral pseudoinverse with cutoff fit_rcond (default 1e-12).
    momentum_fit='legacy' reproduces the C++ iterative fit, with stopping
    tolerance linfit_beta (default 0.1). Supplying linfit_beta selects legacy
    fitting unless momentum_fit is explicitly set. Without either option,
    direct fitting is used. Use the tutorial's linfit_beta=1e-3 to reproduce
    its dynamics: tightening this tolerance changes the propagated basis.
    """
    pars = dyn_params or {}
    q, p, a, s = validate_basis(q, p, a, s)
    ratio, slope, density = momentum_fields(q, p, a, s, coeff_on_surf,
                                                pars.get("node_tol", 1e-14))
    kind = pars.get("mom_calc_type", 1)
    xp = xp_of(q)
    if kind == 0:
        return ratio.imag, ratio.real, xp.zeros_like(q), xp.zeros_like(q)
    if kind == 2:
        return ratio.imag, ratio.real, slope.imag, slope.real
    if kind != 1:
        raise ValueError("mom_calc_type must be 0, 1 or 2")
    weight = density if pars.get("d_weight", 1) else xp.ones_like(density)
    if scalar(weight.max()) == 0:
        return (xp.zeros_like(q),)*4
    weight = weight/weight.max()
    fit_method = pars.get("momentum_fit", "legacy" if "linfit_beta" in pars else "direct")
    if fit_method not in ("direct", "legacy"):
        raise ValueError("momentum_fit must be 'direct' or 'legacy'")
    fitted, slopes = [], []
    for d in range(q.shape[0]):
        design = stack([xp.ones_like(q[d]), q[d]], axis=1)
        normal = design.T @ (weight[:, None]*design)
        rhs = (design.T+0j) @ (weight*ratio[d])
        # normal is real but rhs is complex; promote before PyTorch matmul.
        if fit_method == "legacy":
            fit = _legacy_linear_fit(normal, rhs, pars.get("linfit_beta", .1))
        else:
            fit = metric_solve(normal+0j, rhs, pars.get("fit_rcond", 1e-12))
        fitted.append((design+0j) @ fit)
        slopes.append(xp.ones_like(q[d])*fit[1])
    fitted, slopes = stack(fitted), stack(slopes)
    return fitted.imag, fitted.real, slopes.imag, slopes.real


def time_overlap(nQ, nP, nA, nS, nstate, oQ, oP, oA, oS, ostate):
    """Return <G_new|G_old>, including state masks and QTAG half-widths."""
    return gwp.gwp_overlap_matrix(nQ, nP, nS, nA/2, nstate, oQ, oP, oS, oA/2, ostate)


def propagate_basis(q, p, alp, s, states, coeff, dyn_params, surf_pops):
    """Euler-update the Gaussian basis, then project coefficients onto it.

    q_update_method=1 uses q += dt*invM*momentum; p_update_method=1 replaces
    p with quantum momentum; a_update_method=1 uses alp -= 2*dt*invM*alp*gmom.
    Zero freezes a parameter. Phase evolution is not implemented in legacy
    QTAG: nonzero s_update_method raises instead of silently doing nothing.
    Low-population surfaces synchronize selected parameters to the largest
    population surface (or freeze them); synchronization needs equal counts.
    Stable sorting handles ties. The most populated nonzero state is always
    evolved, even if below decpl_den. A negative new width raises with advice
    to reduce dt; widths and populations are not silently clipped.
    """
    pars = dyn_params
    dt = pars["dt"]
    invM = asarray(pars["iM"], q).reshape(-1, 1)
    target = list(pars.get("states", sorted(set(states))))
    if set(target) != set(states) or len(target) != len(set(target)):
        raise ValueError("states must list each occupied basis label exactly once")
    if len(surf_pops) != len(target):
        raise ValueError("surf_pops must correspond to states")
    if pars.get("s_update_method", 0) != 0:
        raise NotImplementedError("QTAG phase evolution is not implemented; use s_update_method=0")
    updates = [pars.get(k+"_update_method", default) for k, default in (("q", 1), ("p", 1), ("a", 0))]
    sync = [pars.get(k+"_sync_method", default) for k, default in (("q", 1), ("p", 1), ("a", 0), ("s", 0))]
    if any(v not in (0, 1) for v in updates+sync):
        raise ValueError("update and sync methods must be 0 or 1")
    new = [copy(x) for x in (q, p, alp, s)]
    order = sorted(range(len(target)), key=lambda k: -scalar(surf_pops[k]))
    reference = None
    for rank, k in enumerate(order):
        state = target[k]
        idx = indices([i for i, n in enumerate(states) if n == state], q)
        old = [x[:, idx] for x in (q, p, alp, s)]
        pop = scalar(surf_pops[k])
        if pop > pars.get("decpl_den", .25) or (rank == 0 and pop > 0):
            mom, r, gmom, gr = quantum_momentum(*old, coeff.reshape(-1)[idx], pars)
            block = [old[0]+dt*invM*mom if updates[0] else old[0],
                     mom if updates[1] else old[1],
                     old[2]-2*dt*invM*old[2]*gmom if updates[2] else old[2], old[3]]
            if not bool((block[2] > 0).all()):
                raise ValueError("basis update produced nonpositive widths; reduce dt or freeze widths")
            if rank == 0:
                reference = block
        elif reference is not None:
            if any(sync) and old[0].shape != reference[0].shape:
                raise ValueError("synchronization requires equal packet counts on the two states")
            block = [reference[j] if sync[j] else old[j] for j in range(4)]
        else:
            block = old
        for dest, value in zip(new, block):
            dest[:, idx] = value
    nq, np_, na, ns = new
    Snew = time_overlap(nq, np_, na, ns, states, nq, np_, na, ns, states)
    cross = time_overlap(nq, np_, na, ns, states, q, p, alp, s, states)
    cnew = metric_solve(Snew, cross @ coeff.reshape(-1, 1), pars.get("overlap_rcond", 1e-12))
    return nq, np_, na, ns, list(states), cnew
