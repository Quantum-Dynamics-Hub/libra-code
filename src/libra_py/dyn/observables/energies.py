"""
Energy observables computed from TensorStorage.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..hamiltonians import ehrenfest_energy_adi, ehrenfest_energy_dia


def kinetic_energy(storage: Any, traj: Any) -> Any:
    """Return classical kinetic energy for active TBFs."""

    idx = traj.tbf_ids
    p = np.asarray(storage.p[traj.id, idx])
    inv_mass = np.asarray(storage.iM[traj.id, idx])
    return 0.5 * np.sum(p * p * inv_mass, axis=-1)


def active_potential_energy(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """Return active-state potential energy for active TBFs."""

    idx = traj.tbf_ids
    if rep != "adiabatic":
        raise NotImplementedError("active-state energies are implemented for adiabatic rep")

    active = np.asarray(storage.act_states[traj.id, idx], dtype=int)
    energies = np.real(np.diagonal(storage.ham_adi[traj.id, idx], axis1=-2, axis2=-1))
    return energies[np.arange(len(active)), active]


def ehrenfest_potential_energy(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """Return mean-field electronic energy for active TBFs."""

    if rep == "adiabatic":
        return np.real(ehrenfest_energy_adi(storage, traj))
    if rep == "diabatic":
        return np.real(ehrenfest_energy_dia(storage, traj))
    raise ValueError("rep must be 'adiabatic' or 'diabatic'")


def total_energy(
    storage: Any,
    traj: Any,
    rep: str = "adiabatic",
    potential: str = "active",
) -> Any:
    """Return kinetic plus selected potential energy for active TBFs."""

    if potential == "active":
        pe = active_potential_energy(storage, traj, rep)
    elif potential == "ehrenfest":
        pe = ehrenfest_potential_energy(storage, traj, rep)
    else:
        raise ValueError("potential must be 'active' or 'ehrenfest'")
    return kinetic_energy(storage, traj) + pe
