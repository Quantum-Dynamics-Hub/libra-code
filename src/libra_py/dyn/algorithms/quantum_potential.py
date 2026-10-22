# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Quantum momentum and logarithmic-gradient estimation for Gaussian expansions.

Stateless NumPy/PyTorch kernels accept (ndof, npackets) arrays. Direct spectral
and tutorial-compatible iterative fits share the Gaussian wavefunction kernels.
These routines estimate fields and slopes; basis propagation and coefficient
projection remain in propagation.gaussian. No scalar quantum-potential energy
is computed here.
"""
import numpy as np
from ..gwp.packets import validate_basis, momentum_fields
from ..utils.linalg import metric_solve
from ..backends.arrays import xp_of, stack, scalar


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
    """Return momentum, real log-gradient, and their slopes (ndof,npackets).

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


