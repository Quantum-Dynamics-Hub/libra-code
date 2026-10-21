# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Initial Gaussian bases and wavefunction projection for array-based QTAG."""
import numpy as np
from ...gwp import integrals as gwp
from . import core
from ._backend import asarray, xp_of, scalar


def _convert(values, params):
    backend = params.get("backend", "numpy")
    if backend == "numpy":
        return np.asarray(values, dtype=params.get("dtype", "float64"))
    if backend != "torch":
        raise ValueError("backend must be numpy or torch")
    import torch
    dtype = params.get("dtype", torch.float64)
    if isinstance(dtype, str):
        dtype = getattr(torch, dtype)
    return torch.as_tensor(values, dtype=dtype, device=params.get("device", "cpu"))


def _parameters(params):
    q0 = np.asarray(params.get("wfc_q0", [0.]), dtype=float)
    if q0.ndim != 1 or len(q0) == 0:
        raise ValueError("wfc_q0 must be a nonempty vector")
    d = len(q0)
    def vector(key, default):
        x = np.broadcast_to(np.asarray(params.get(key, default), dtype=float), (d,)).copy()
        if not np.isfinite(x).all():
            raise ValueError(key+" must be finite")
        return x
    p0, a0, s0 = vector("wfc_p0", 0.), vector("wfc_a0", 1.), vector("wfc_s0", 0.)
    scale = vector("alp_scl", 8.)
    if np.any(a0 <= 0) or np.any(scale <= 0):
        raise ValueError("initial widths and alp_scl must be positive")
    states = list(params.get("states", [0, 1]))
    core._states(states, len(states))
    if len(set(states)) != len(states):
        raise ValueError("states must be unique")
    return q0, p0, a0, s0, scale, states


def _assemble(positions, momenta, a0, s0, scale, states, params):
    n = positions.shape[1]
    # All states start from the same basis, giving a clear synchronization map.
    q = np.tile(positions, (1, len(states)))
    p = np.tile(momenta, (1, len(states)))
    a = np.broadcast_to((a0*scale)[:, None], q.shape).copy()
    s = np.broadcast_to(s0[:, None], q.shape).copy()
    return (q.shape[1], *[_convert(x, params) for x in (q, p, a, s)],
            [state for state in states for _ in range(n)])


def grid(params):
    """Place a Cartesian basis on each state; return ntraj,q,p,alp,s,labels.

    grid_dims gives packet counts per DOF (default 5 each). For 0<rho_cut<1,
    bounds are q0 +/- sqrt(-log(rho_cut)/(2*wfc_a0)). Otherwise explicit
    grid_min/max are used. Basis QTAG widths are wfc_a0*alp_scl, as in the
    legacy tutorial. backend='torch' optionally selects device/dtype.
    """
    q0, p0, a0, s0, scale, states = _parameters(params)
    dims = params.get("grid_dims", [5]*len(q0))
    if len(dims) != len(q0) or any(int(n) != n or n < 1 for n in dims):
        raise ValueError("grid_dims needs one positive integer per DOF")
    rho = params.get("rho_cut", 1e-12)
    if rho > 0:
        if rho >= 1:
            raise ValueError("positive rho_cut must be less than one")
        delta = np.sqrt(-np.log(rho)/(2*a0))
        lo, hi = q0-delta, q0+delta
    else:
        lo = np.broadcast_to(params.get("grid_min", -1.), q0.shape)
        hi = np.broadcast_to(params.get("grid_max", 1.), q0.shape)
        if np.any(hi <= lo):
            raise ValueError("grid_max must exceed grid_min")
    axes = [np.linspace(l, h, int(n)) if n > 1 else np.array([(l+h)/2])
            for l, h, n in zip(lo, hi, dims)]
    q = np.stack(np.meshgrid(*axes, indexing="ij")).reshape(len(q0), -1)
    return _assemble(q, np.broadcast_to(p0[:, None], q.shape), a0, s0, scale, states, params)


def gaussian(params):
    """Sample a reproducible Gaussian basis, cloned onto each electronic state.

    prod(grid_dims) is the number per state. rho_cut is the sampling standard
    deviation (legacy meaning), default 1.0 here; seed controls a local RNG.
    Coordinates and momenta are sorted per DOF as in the legacy initializer.
    Random draws use NumPy before conversion to the requested backend.
    """
    q0, p0, a0, s0, scale, states = _parameters(params)
    dims = params.get("grid_dims", [5])
    if not dims or any(int(n) != n or n < 1 for n in dims):
        raise ValueError("grid_dims must contain positive integers")
    n = int(np.prod(dims))
    sigma = params.get("rho_cut", 1.)
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("Gaussian sampling rho_cut must be positive")
    rng = np.random.default_rng(params.get("seed", 0))
    q = np.sort(rng.normal(q0[:, None], sigma, (len(q0), n)), axis=1)
    p = np.sort(rng.normal(p0[:, None], sigma, (len(q0), n)), axis=1)
    return _assemble(q, p, a0, s0, scale, states, params)


def gbfs(params):
    """Select grid (init_placement=0) or Gaussian sampling (1)."""
    method = params.get("init_placement", 0)
    if method not in (0, 1):
        raise ValueError("init_placement must be 0 or 1")
    return (grid if method == 0 else gaussian)(params)


def wfc_centers(params):
    """Return q0,p0,a0,s0 as (ndof,1) arrays in the chosen backend.

    IMPORTANT: wfc_a0 uses the original GWP exponent -a0*(r-q0)^2, while
    basis widths use the QTAG exponent -alp*(r-q)^2/2. This distinction
    preserves the initial wavepacket from the original tutorial.
    """
    q, p, a, s, _, _ = _parameters(params)
    return tuple(_convert(x[:, None], params) for x in (q, p, a, s))


def coeffs(aQ, aP, aA, aS, astate, bQ, bP, bA, bS, bstate,
           rcond=1e-12, normalize=True):
    """Project one initial GWP onto a state-labeled QTAG basis.

    Returns (ntraj,1) complex coefficients. aA uses GWP widths; bA uses
    QTAG widths. Solve in the retained positive overlap eigenspace instead
    of explicitly inverting S. normalize=True normalizes the projected
    state; False reproduces the legacy unnormalized projection.
    """
    aQ, aP, aA, aS = [asarray(x, bQ) for x in (aQ, aP, aA, aS)]
    if aQ.shape != (bQ.shape[0], 1):
        raise ValueError("initial wavefunction must be one packet")
    rhs = gwp.gwp_overlap_matrix(bQ, bP, bS, bA/2, bstate, aQ, aP, aS, aA, astate)
    S = gwp.gwp_overlap_matrix(bQ, bP, bS, bA/2, bstate, bQ, bP, bS, bA/2, bstate)
    c = core.metric_solve(S, rhs, rcond)
    if normalize:
        norm = (c.conj().T @ S @ c)[0, 0].real
        if scalar(norm) <= 0:
            raise ValueError("initial packet has no projection onto the selected states")
        c = c/xp_of(c).sqrt(norm)
    return c
