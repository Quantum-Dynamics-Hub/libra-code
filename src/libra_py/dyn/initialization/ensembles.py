# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Shared nuclear sampling and initialization of independent/coherent TBF ensembles.

Independent paths normalize electronic amplitudes per path. Coherent Gaussian
expansions instead project and normalize the whole wavefunction with its overlap
metric. Both use the same nuclear sampling and storage assignment helpers.
"""
from __future__ import annotations

import numpy as np

from ..backends import backend as default_backend
from ..core.storage import TensorStorage
from ..core.trajectory import Trajectory
from ..gwp import integrals as gwp
from ..gwp.packets import validate_states
from ..utils.linalg import metric_solve
from ..backends.arrays import asarray, xp_of, scalar


def sample_gaussian_initial_conditions(
    ntraj: int,
    q_mean,
    q_sigma,
    p_mean,
    p_sigma,
    ndof: int = 1,
    rng=None,
    *,
    draw_order="trajectory",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample independent Gaussian nuclear coordinates and momenta.

    The returned arrays have shape ``(ntraj, ndof)`` and are suitable for
    assigning to ``TensorStorage.q[traj.id, traj.tbf_ids]`` and ``p``. Scalar
    means/widths are broadcast over degrees of freedom.

    ``draw_order="trajectory"`` retains the original RNG draw order and seeded
    results. ``"dof"`` draws each coordinate's samples consecutively, matching
    the historical QTAG sampler. Both return (ntraj,ndof); neither sorts samples.
    Coordinates are drawn before momenta in both modes.
    """

    rng = rng or np.random.default_rng()
    q_mean = np.broadcast_to(np.asarray(q_mean, dtype=float), (ndof,))
    q_sigma = np.broadcast_to(np.asarray(q_sigma, dtype=float), (ndof,))
    p_mean = np.broadcast_to(np.asarray(p_mean, dtype=float), (ndof,))
    p_sigma = np.broadcast_to(np.asarray(p_sigma, dtype=float), (ndof,))
    if draw_order == "trajectory":
        q = rng.normal(q_mean, q_sigma, size=(ntraj, ndof))
        p = rng.normal(p_mean, p_sigma, size=(ntraj, ndof))
    elif draw_order == "dof":
        q = rng.normal(q_mean[:, None], q_sigma[:, None], size=(ndof, ntraj)).T
        p = rng.normal(p_mean[:, None], p_sigma[:, None], size=(ndof, ntraj)).T
    else:
        raise ValueError("draw_order must be 'trajectory' or 'dof'")
    return q, p


def make_independent_trajectory_ensemble(
    q,
    p,
    masses,
    amplitudes,
    backend=default_backend,
    traj_id: int = 0,
    active_state: int = 0,
    der_lvl: int = 1,
) -> tuple[TensorStorage, Trajectory]:
    """
    Build a storage/trajectory pair for independent mixed quantum-classical paths.

    ``q`` and ``p`` have shape ``(ntraj, ndof)``. ``masses`` may be scalar or
    length ``ndof``. ``amplitudes`` may be a single state vector of shape
    ``(nstates,)`` or one vector per path with shape ``(ntraj, nstates)``.
    """

    q = _as_2d(q, "q")
    p = _as_2d(p, "p")
    if q.shape != p.shape:
        raise ValueError("q and p must have the same shape")
    npaths, ndof = q.shape

    amplitudes = np.asarray(amplitudes, dtype=complex)
    if amplitudes.ndim == 1:
        amplitudes = np.broadcast_to(amplitudes, (npaths, amplitudes.shape[0])).copy()
    if amplitudes.shape[0] != npaths:
        raise ValueError("amplitudes must have one row per trajectory path")
    nstates = amplitudes.shape[-1]

    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=ndof,
        nstates=nstates,
        ntbf_initial=npaths,
        ntbf_capacity=npaths,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=der_lvl)

    traj = Trajectory(traj_id)
    traj.tbf_ids = list(range(npaths))
    idx = traj.tbf_ids
    _initialize_nuclear(storage, traj, idx, q, p,
                        1.0 / np.broadcast_to(np.asarray(masses, dtype=float), (ndof,)))
    storage.ampl_adi[traj.id, idx] = _normalized_amplitudes(amplitudes)
    storage.act_states[traj.id, idx] = int(active_state)
    return storage, traj


def _as_2d(value, name: str):
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2:
        raise ValueError(f"{name} must have shape (ntraj, ndof)")
    return arr


def _normalized_amplitudes(amplitudes):
    norm = np.linalg.norm(amplitudes, axis=-1)
    if np.any(norm == 0.0):
        raise ValueError("amplitudes must be nonzero")
    return amplitudes / norm[:, None]


def _initialize_nuclear(storage, traj, slots, q, p, inverse_masses):
    """Populate common nuclear fields for either ensemble representation.

    Inputs use (npackets,ndof). Conversion follows the storage dtype/device;
    inverse masses must be common, finite and positive for each DOF.
    Electronic normalization belongs to the respective ensemble initializer.
    """
    like = storage.q
    q, p = asarray(q, like), asarray(p, like)
    shape = (len(slots), storage.ndof)
    if q.shape != shape or p.shape != shape:
        raise ValueError(f"q and p must have shape {shape}")
    mass = asarray(inverse_masses, like).reshape(-1)
    if mass.shape != (storage.ndof,) or not bool(xp_of(mass).isfinite(mass).all()) or not bool((mass > 0).all()):
        raise ValueError('one finite positive inverse mass per DOF is required')
    storage.q[traj.id, slots] = q
    storage.p[traj.id, slots] = p
    storage.iM[traj.id, slots] = mass


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
    validate_states(states, len(states))
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
    q, p = sample_gaussian_initial_conditions(
        n, q0, sigma, p0, sigma, ndof=len(q0), rng=rng, draw_order="dof")
    q, p = np.sort(q, axis=0).T, np.sort(p, axis=0).T
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
    c = metric_solve(S, rhs, rcond)
    if normalize:
        norm = (c.conj().T @ S @ c)[0, 0].real
        if scalar(norm) <= 0:
            raise ValueError("initial packet has no projection onto the selected states")
        c = c/xp_of(c).sqrt(norm)
    return c


def make_gaussian_ensemble(params, inverse_masses, nstates=None):
    """Initialize and register one or several coherent Gaussian wavefunctions.

    Parameters
    ----------
    params : dict or sequence of dict
        Each dictionary uses grid/gbfs options above, plus ``initial_state``
        (default first state label), ``normalize`` (True), and ``overlap_rcond``
        (1e-12). A sequence creates independent Trajectories, each owning its
        own coherent set of TBFs; packet counts may differ. All must share
        ndof, backend, dtype and device. Sampling uses a local NumPy RNG.
    inverse_masses : array-like, (ndof,)
        Common positive inverse nuclear masses, in atomic units.
    nstates : int, optional
        Electronic storage dimension; defaults to largest label plus one.

    Returns
    -------
    TrajectoryManager
        Owns TensorStorage and all registered Trajectory/TBF objects. Live
        q,p,alpha,phase use (ntraj,capacity,ndof), amplitudes use
        (ntraj,capacity,nstates). Only the labeled amplitude is nonzero.
        No per-packet normalization is applied to the projected coefficients.
    """
    from ..backends import NumpyBackend
    from ..backends.torch_backend import TorchBackend
    from ..core.manager import TrajectoryManager
    from ..gwp.storage import write_coefficients

    specifications = [dict(params)] if isinstance(params, dict) else [dict(p) for p in params]
    if not specifications:
        raise ValueError('at least one trajectory specification is required')
    bases = [gbfs(p) for p in specifications]
    ndof = bases[0][1].shape[0]
    first = specifications[0]
    backend = first.get('backend', 'numpy')
    if any(b[1].shape[0] != ndof or p.get('backend', 'numpy') != backend
           or str(p.get('dtype', 'float64')) != str(first.get('dtype', 'float64'))
           or p.get('device', 'cpu') != first.get('device', 'cpu')
           for p, b in zip(specifications, bases)):
        raise ValueError('ensemble members must share ndof, backend, dtype and device')
    count = max(max(b[-1]) for b in bases)+1
    nstates = count if nstates is None else nstates
    if int(nstates) != nstates or nstates < count:
        raise ValueError('nstates must include every packet label')
    bd = (TorchBackend(first.get('device', 'cpu'), first.get('dtype', 'float64'))
          if backend == 'torch' else NumpyBackend(first.get('dtype', 'float64')))
    capacity = max(b[0] for b in bases)
    storage = TensorStorage(bd, len(bases), ndof, int(nstates), capacity, capacity)
    storage.allocate_gaussian_vars()
    storage.alive[:] = False
    manager = TrajectoryManager(storage)
    for spec, (n, q, p, a, s, labels) in zip(specifications, bases):
        traj = manager.create_trajectory(metadata={'representation': 'state-labeled Gaussians'})
        for slot in range(n):
            manager.create_tbf(traj, slot)
        _initialize_nuclear(storage, traj, list(range(n)), q.T, p.T, inverse_masses)
        storage.gaussian_alpha[traj.id, :n] = a.T
        storage.gaussian_phase[traj.id, :n] = s.T
        storage.act_states_dia[traj.id, :n] = bd.asarray(labels, dtype=int)
        initial = int(spec.get('initial_state', labels[0]))
        q0, p0, a0, s0 = wfc_centers(spec)
        c = coeffs(q0, p0, a0, s0, [initial], q, p, a, s, labels,
                   spec.get('overlap_rcond', 1e-12), spec.get('normalize', True))
        write_coefficients(storage, traj, list(range(n)), labels, c)
    return manager
