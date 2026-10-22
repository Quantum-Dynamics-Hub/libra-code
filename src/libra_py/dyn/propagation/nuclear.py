"""
Classical nuclear propagation kernels.

The functions here operate on active TensorStorage slices and intentionally
avoid Hamiltonian/electronic decisions.
"""

from __future__ import annotations

from typing import Any


def drift_coordinates(q, momentum, inverse_mass, dt):
    """Return q + dt * inverse_mass * momentum without modifying inputs.

    Arrays must have broadcast-compatible layouts and use the same backend.
    Supports NumPy/PyTorch and arbitrary batch/DOF axes. Momentum may be a
    classical stored momentum or a quantum-momentum field supplied by QTAG.
    """
    return q + dt * inverse_mass * momentum


def velocity(storage: Any, traj: Any) -> Any:
    """Return active nuclear velocities v = p / M = p * iM."""

    idx = traj.tbf_ids
    return storage.p[traj.id, idx] * storage.iM[traj.id, idx]


def drift(storage: Any, traj: Any, dt: float) -> Any:
    """Advance active nuclear coordinates by one drift substep."""

    idx = traj.tbf_ids
    storage.q[traj.id, idx] = drift_coordinates(
        storage.q[traj.id, idx], storage.p[traj.id, idx],
        storage.iM[traj.id, idx], dt,
    )
    return storage.q[traj.id, idx]


def kick(storage: Any, traj: Any, dt: float, forces: Any = None) -> Any:
    """Advance active nuclear momenta by one force kick substep."""

    idx = traj.tbf_ids
    if forces is None:
        forces = storage.f[traj.id, idx]
    storage.p[traj.id, idx] = storage.p[traj.id, idx] + dt * forces
    return storage.p[traj.id, idx]


def velocity_verlet_step(
    storage: Any,
    traj: Any,
    dt: float,
    force_fn,
    before_force_fn=None,
) -> Any:
    """
    Full velocity-Verlet step.

    `force_fn` must write/return forces at the current coordinates. The
    optional `before_force_fn` hook is called after the drift and before the
    second force evaluation, which is where Hamiltonians are normally updated.
    """

    kick(storage, traj, 0.5 * dt)
    drift(storage, traj, dt)
    if before_force_fn is not None:
        before_force_fn()
    forces = force_fn()
    kick(storage, traj, 0.5 * dt, forces)
    return storage
