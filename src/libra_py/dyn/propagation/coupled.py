"""
Coupled electron-nuclear propagation helpers.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..hamiltonians import ehrenfest_forces_adi, ehrenfest_forces_dia


def density_matrix(amplitudes: Any) -> Any:
    """Return |C><C| for an active trajectory/TBF amplitude batch."""

    c = np.asarray(amplitudes)
    return np.einsum("...i,...j->...ij", c, np.conjugate(c))


def update_density(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """Update stored density matrices from the current amplitudes."""

    idx = traj.tbf_ids
    if rep == "adiabatic":
        storage.dm_adi[traj.id, idx] = density_matrix(storage.ampl_adi[traj.id, idx])
        return storage.dm_adi[traj.id, idx]
    if rep == "diabatic":
        storage.dm_dia[traj.id, idx] = density_matrix(storage.ampl_dia[traj.id, idx])
        return storage.dm_dia[traj.id, idx]
    raise ValueError("rep must be 'adiabatic' or 'diabatic'")


def state_specific_forces(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """
    Return active state-specific forces.

    Currently implemented for adiabatic gradients, where
    F_I = -d E_I / dR.
    """

    idx = traj.tbf_ids
    if rep != "adiabatic":
        raise NotImplementedError("state-specific forces are implemented for adiabatic rep")
    if storage.d1ham_adi is None:
        raise AttributeError("TensorStorage field 'd1ham_adi' has not been allocated")

    active = np.asarray(storage.act_states[traj.id, idx], dtype=int)
    d1 = np.asarray(storage.d1ham_adi[traj.id, idx])
    gradients = np.diagonal(d1, axis1=-2, axis2=-1)
    rows = np.arange(len(active))[:, None]
    dofs = np.arange(gradients.shape[1])[None, :]
    return -np.real(gradients[rows, dofs, active[:, None]])


def ehrenfest_forces(
    storage: Any,
    traj: Any,
    rep: str = "adiabatic",
    option: int = 0,
    transform: Any = None,
    gamma: float | Any = 0.0,
) -> Any:
    """Return active Ehrenfest mean-field forces."""

    if rep == "adiabatic":
        return np.real(ehrenfest_forces_adi(
            storage,
            traj,
            option=option,
            transform=transform,
            gamma=gamma,
        ))
    if rep == "diabatic":
        return np.real(ehrenfest_forces_dia(storage, traj, option=option))
    raise ValueError("rep must be 'adiabatic' or 'diabatic'")
