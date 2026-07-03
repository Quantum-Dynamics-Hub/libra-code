"""
Population observables computed from TensorStorage.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def amplitude_populations(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """Return per-state populations |C_i|^2 for active TBFs."""

    amplitudes = _amplitudes(storage, traj, rep)
    return np.real(np.conjugate(amplitudes) * amplitudes)


def density_populations(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """Return per-state populations from the diagonal of stored density matrices."""

    idx = traj.tbf_ids
    if rep == "adiabatic":
        density = storage.dm_adi[traj.id, idx]
    elif rep == "diabatic":
        density = storage.dm_dia[traj.id, idx]
    else:
        raise ValueError("rep must be 'adiabatic' or 'diabatic'")
    return np.real(np.diagonal(density, axis1=-2, axis2=-1))


def active_state_counts(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """Count active states over active TBFs for one trajectory."""

    idx = traj.tbf_ids
    states = (
        storage.act_states[traj.id, idx]
        if rep == "adiabatic"
        else storage.act_states_dia[traj.id, idx]
    )
    return np.bincount(np.asarray(states, dtype=int), minlength=storage.nstates)


def mean_populations(storage: Any, traj: Any, rep: str = "adiabatic") -> Any:
    """Return TBF-averaged amplitude populations for one trajectory."""

    pops = amplitude_populations(storage, traj, rep)
    if pops.shape[0] == 0:
        return np.zeros(storage.nstates)
    return np.mean(pops, axis=0)


def _amplitudes(storage: Any, traj: Any, rep: str) -> Any:
    idx = traj.tbf_ids
    if rep == "adiabatic":
        return storage.ampl_adi[traj.id, idx]
    if rep == "diabatic":
        return storage.ampl_dia[traj.id, idx]
    raise ValueError("rep must be 'adiabatic' or 'diabatic'")
