from __future__ import annotations

import numpy as np

from ..backends import backend as default_backend
from ..core.storage import TensorStorage
from ..core.trajectory import Trajectory


def sample_gaussian_initial_conditions(
    ntraj: int,
    q_mean,
    q_sigma,
    p_mean,
    p_sigma,
    ndof: int = 1,
    rng=None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample independent Gaussian nuclear coordinates and momenta.

    The returned arrays have shape ``(ntraj, ndof)`` and are suitable for
    assigning to ``TensorStorage.q[traj.id, traj.tbf_ids]`` and ``p``. Scalar
    means/widths are broadcast over degrees of freedom.
    """

    rng = rng or np.random.default_rng()
    q_mean = np.broadcast_to(np.asarray(q_mean, dtype=float), (ndof,))
    q_sigma = np.broadcast_to(np.asarray(q_sigma, dtype=float), (ndof,))
    p_mean = np.broadcast_to(np.asarray(p_mean, dtype=float), (ndof,))
    p_sigma = np.broadcast_to(np.asarray(p_sigma, dtype=float), (ndof,))
    q = rng.normal(q_mean, q_sigma, size=(ntraj, ndof))
    p = rng.normal(p_mean, p_sigma, size=(ntraj, ndof))
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
    storage.q[traj.id, idx] = q
    storage.p[traj.id, idx] = p
    storage.iM[traj.id, idx] = 1.0 / np.broadcast_to(np.asarray(masses, dtype=float), (ndof,))
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
