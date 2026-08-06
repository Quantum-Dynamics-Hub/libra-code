"""
Electronic propagation module.

Responsibilities:
------------------
- TDSE integrators (pure kernels)
- Basis transformations (local diabatization)
- TDSE orchestration (given TensorStorage)
- Backend-agnostic linear algebra

NOT responsible for:
--------------------
- Hamiltonian construction
- NAC computation
- Force evaluation
- Trajectory management
"""

from __future__ import annotations

from typing import Callable, Optional, Literal

from ..backends import backend as backend_default


# ============================================================
# Types
# ============================================================

Backend = backend_default.__class__

Representation = Literal["adiabatic", "diabatic"]
HamiltonianType = Literal["hamiltonian", "vibronic"]
Propagator = Callable


# ============================================================
# Core propagators
# ============================================================

def _apply_matrix_to_state(A, C, backend=backend_default):
    if getattr(C, "ndim", 0) == getattr(A, "ndim", 0) - 1:
        return backend.einsum("...ij,...j->...i", A, C)
    return backend.matmul(A, C)


def exp_propagator(C, H, dt, backend=backend_default):
    """
    Standard unitary TDSE propagation:
        C(t+dt) = exp(-i H dt) C(t)
    """
    U = backend.expm(-1j * H * dt)
    return _apply_matrix_to_state(U, C, backend)


def euler_propagator(C, H, dt, backend=backend_default):
    """
    First-order Euler TDSE integrator (debug / weak coupling only)
    """
    HC = _apply_matrix_to_state(H, C, backend)
    return C - 1j * dt * HC


def split_step_propagator(C, H_prev, H, dt, backend=backend_default):
    """
    Symmetric 2-point propagator:
        exp(-i H_prev dt/2) exp(-i H dt/2)
    """
    U1 = backend.expm(-1j * H_prev * dt * 0.5)
    U2 = backend.expm(-1j * H * dt * 0.5)

    return _apply_matrix_to_state(
        U2,
        _apply_matrix_to_state(U1, C, backend),
        backend
    )


# ============================================================
# Basis transformations
# ============================================================

def apply_local_diabatization(C, H, T, backend=backend_default):
    """
    Local diabatization transform:

        H' = T† H T
        C' = T† C
    """
    Tdag = backend.conjugate_transpose(T)

    H_rot = backend.matmul(Tdag, backend.matmul(H, T))
    C_rot = _apply_matrix_to_state(Tdag, C, backend)

    return C_rot, H_rot


# ============================================================
# TDSE orchestration
# ============================================================

def tdse_step(
    traj,
    storage,
    state_or_dt=None,
    dt=None,
    backend=backend_default,
    propagator: Propagator = exp_propagator,
    rep: Representation = "adiabatic",
    hamiltonian_type: HamiltonianType = "vibronic",
    T: Optional[object] = None,
    previous_state: Optional[object] = None,
):
    """
    Full electronic propagation step.

    Parameters
    ----------
    traj :
        Trajectory object (defines active TBF indices)

    storage :
        TensorStorage (C, R, etc.)

    state_or_dt :
        Either the timestep dt, or a deprecated Hamiltonian snapshot when dt is
        supplied separately.

    dt :
        Time step. If omitted, state_or_dt is treated as dt.

    backend :
        linear algebra backend (NumPy / PyTorch / JAX)

    propagator :
        TDSE integrator

    rep :
        "adiabatic" or "diabatic"

    hamiltonian_type :
        "hamiltonian" or "vibronic"

    T :
        optional local diabatization transform

    previous_state :
        optional previous Hamiltonian snapshot used by two-point propagators.
    """

    state = None
    if dt is None:
        dt = state_or_dt
    else:
        state = state_or_dt

    idx = traj.tbf_ids

    # --------------------------------------------------------
    # 1. Load electronic amplitudes
    # --------------------------------------------------------
    if rep == "adiabatic":
        C = storage.ampl_adi[traj.id, idx]
    else:
        C = storage.ampl_dia[traj.id, idx]

    # --------------------------------------------------------
    # 2. Select Hamiltonian representation
    # --------------------------------------------------------
    H = _hamiltonian_slice(
        storage,
        traj,
        rep,
        hamiltonian_type,
        state=state,
    )
    H_prev = _previous_hamiltonian_slice(
        previous_state if previous_state is not None else state,
        rep,
    )

    # --------------------------------------------------------
    # 4. Local diabatization (basis transform)
    # --------------------------------------------------------
    if T is not None:
        C, H = apply_local_diabatization(C, H, T, backend)

    # --------------------------------------------------------
    # 5. Propagation
    # --------------------------------------------------------
    if propagator is split_step_propagator and H_prev is not None:
        C_new = propagator(C, H_prev, H, dt, backend)
    else:
        C_new = propagator(C, H, dt, backend)

    # --------------------------------------------------------
    # 6. Write back
    # --------------------------------------------------------
    if rep == "adiabatic":
        storage.ampl_adi[traj.id, idx] = C_new
    else:
        storage.ampl_dia[traj.id, idx] = C_new

    return C_new


# ============================================================
# Convenience wrapper
# ============================================================

def propagate_electronic(*args, **kwargs):
    """
    Alias for tdse_step for backward compatibility.
    """
    return tdse_step(*args, **kwargs)


def _hamiltonian_slice(
    storage,
    traj,
    rep: Representation,
    hamiltonian_type: HamiltonianType,
    state=None,
):
    """Return the Hamiltonian matrix batch for one trajectory."""

    if state is not None and state is not storage:
        return _legacy_state_matrix(state, rep, hamiltonian_type)

    if rep == "adiabatic":
        field = "hvib_adi" if hamiltonian_type == "vibronic" else "ham_adi"
    else:
        field = "hvib_dia" if hamiltonian_type == "vibronic" else "ham_dia"

    matrix = getattr(storage, field)
    if matrix is None:
        raise AttributeError(
            f"TensorStorage field '{field}' has not been allocated"
        )
    return matrix[traj.id, traj.tbf_ids]


def _legacy_state_matrix(
    state,
    rep: Representation,
    hamiltonian_type: HamiltonianType,
):
    """Read from the former HamiltonianState shape when callers still pass it."""

    if isinstance(state, dict):
        suffix = "adi" if rep == "adiabatic" else "dia"
        prefix = "hvib" if hamiltonian_type == "vibronic" else "ham"
        return state[f"{prefix}_{suffix}"]

    if rep == "adiabatic":
        return state.Hvib_adi if hamiltonian_type == "vibronic" else state.H_adi
    return state.Hvib_dia if hamiltonian_type == "vibronic" else state.H_dia


def _previous_hamiltonian_slice(state, rep: Representation):
    if state is None:
        return None
    if isinstance(state, dict):
        return state[f"hvib_{'adi' if rep == 'adiabatic' else 'dia'}"]
    if rep == "adiabatic":
        return getattr(state, "H_adi_prev", None)
    return getattr(state, "H_dia_prev", None)
