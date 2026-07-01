"""
Electronic propagation module.

Responsibilities:
------------------
- TDSE integrators (pure kernels)
- Basis transformations (local diabatization)
- TDSE orchestration (given HamiltonianState)
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
    state,
    dt,
    backend=backend_default,
    propagator: Propagator = exp_propagator,
    rep: Representation = "adiabatic",
    hamiltonian_type: HamiltonianType = "vibronic",
    T: Optional[object] = None,
):
    """
    Full electronic propagation step.

    Parameters
    ----------
    traj :
        Trajectory object (defines active TBF indices)

    storage :
        TensorStorage (C, R, etc.)

    state :
        HamiltonianState from HamiltonianEngine.evaluate()

    dt :
        time step

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
    """

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
    if rep == "adiabatic":
        H = state.H_adi
        H_prev = getattr(state, "H_adi_prev", None)
    else:
        H = state.H_dia
        H_prev = getattr(state, "H_dia_prev", None)

    # --------------------------------------------------------
    # 3. Optional vibronic Hamiltonian
    # --------------------------------------------------------
    if hamiltonian_type == "vibronic":
        if rep == "adiabatic":
            H = state.Hvib_adi
        else:
            H = state.Hvib_dia

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
