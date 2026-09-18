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

import numpy as np

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


def _propagate_fixed_hamiltonian(
    C,
    H,
    dt,
    nsubsteps=1,
    backend=backend_default,
    propagator: Propagator = exp_propagator,
):
    """Propagate electronic states while holding the Hamiltonian fixed.

    Parameters
    ----------
    C : backend array, shape ``(..., nstates)`` or ``(..., nstates, nvec)``
        Electronic state vector(s) at the beginning of the interval. Leading
        dimensions identify independent trajectories or other batches. The
        penultimate dimension is the electronic-state dimension when multiple
        vectors are stored in the final dimension.
    H : backend array, shape ``(..., nstates, nstates)``
        Fixed electronic Hamiltonian. Its leading dimensions must be
        broadcast-compatible with those of ``C``.
    dt : float
        Total propagation interval, in the time units reciprocal to the
        energy units of ``H`` (atomic units in dynamics calculations).
    nsubsteps : int, optional
        Number of equal intervals into which ``dt`` is divided. Expected to be
        a positive integer; validation is performed by the public caller.
    backend : linear-algebra backend, optional
        Backend object supplying the array operations needed by
        ``propagator``. The default is the configured Libra backend.
    propagator : callable, optional
        One-step function with signature ``(C, H, dt, backend) -> C_new``.

    Returns
    -------
    backend array, same shape as ``C``
        Electronic state vector(s) after ``nsubsteps`` applications over the
        complete interval ``dt``. The input arrays are not explicitly
        modified by this helper.
    """

    state = C
    substep_dt = dt / nsubsteps
    for _ in range(nsubsteps):
        state = propagator(state, H, substep_dt, backend)
    return state


def _rotate_matrix(matrix, transform, transform_h, backend=backend_default):
    """Express a matrix in the old dynamically tracked adiabatic basis.

    Parameters
    ----------
    matrix : backend array, shape ``(..., nstates, nstates)``
        Matrix represented in the current raw adiabatic basis. In this module
        it is normally an electronic or vibronic Hamiltonian.
    transform : backend array, shape ``(..., nstates, nstates)``
        Old-to-current basis transformation ``T``. Its columns express old
        dynamically tracked state labels in the current raw basis.
    transform_h : backend array, shape ``(..., nstates, nstates)``
        Conjugate transpose ``T†`` of ``transform``. It is accepted explicitly
        so callers that rotate several matrices can compute it once.
    backend : linear-algebra backend, optional
        Backend object providing batched ``matmul``.

    Returns
    -------
    backend array, shape ``(..., nstates, nstates)``
        The matrix ``T† @ matrix @ T`` in the old tracked basis. Its dtype is
        determined by the input arrays and backend, and is normally complex.
    """

    return backend.matmul(transform_h, backend.matmul(matrix, transform))


def interpolated_hamiltonian_propagator(
    C,
    H_previous,
    H_current,
    dt,
    nsubsteps=1,
    backend=backend_default,
    propagator: Propagator = exp_propagator,
):
    """Propagate through a Hamiltonian that varies linearly in time.

    Each substep uses the Hamiltonian at its temporal midpoint,

    ``H(alpha) = (1 - alpha) H_previous + alpha H_current``.

    Both endpoint matrices must be expressed in the same basis. The midpoint
    rule is second-order accurate and is exact for scalar/commuting endpoint
    Hamiltonians.

    Parameters
    ----------
    C : backend array, shape ``(..., nstates)`` or ``(..., nstates, nvec)``
        Electronic state vector(s) at the previous endpoint. Leading
        dimensions may batch independent trajectories.
    H_previous : backend array, shape ``(..., nstates, nstates)``
        Hamiltonian at the beginning of the interval. For adiabatic dynamics,
        it contains the previous-endpoint adiabatic energies and any coupling
        terms required by the selected electronic equation of motion.
    H_current : backend array, shape ``(..., nstates, nstates)``
        Hamiltonian at the end of the interval, expressed in the same basis
        and state ordering as ``H_previous``. Its entries are linearly
        interpolated with those of ``H_previous``.
    dt : float
        Total endpoint-to-endpoint propagation interval, in the time units
        reciprocal to the Hamiltonian energy units.
    nsubsteps : int, optional
        Positive number of equal electronic substeps. Substep ``k`` uses the
        midpoint fraction ``alpha = (k + 1/2) / nsubsteps``.
    backend : linear-algebra backend, optional
        Backend object used for array and matrix operations.
    propagator : callable, optional
        One-step function with signature ``(C, H, dt, backend) -> C_new``.
        The default applies the matrix exponential ``exp(-i H dt)``.

    Returns
    -------
    backend array, same shape as ``C``
        Electronic state vector(s) at the current endpoint. The returned dtype
        is determined by ``C``, the Hamiltonians, and the selected propagator;
        TD-SE propagation normally returns a complex array. Endpoint
        Hamiltonians and ``C`` are not explicitly modified.

    Raises
    ------
    ValueError
        If ``nsubsteps`` is less than one.
    """

    nsubsteps = int(nsubsteps)
    if nsubsteps < 1:
        raise ValueError("nsubsteps must be positive")
    state = C
    substep_dt = dt / nsubsteps
    for substep in range(nsubsteps):
        alpha = (substep + 0.5) / nsubsteps
        hamiltonian = (1.0 - alpha) * H_previous + alpha * H_current
        state = propagator(state, hamiltonian, substep_dt, backend)
    return state


def propagate_electronic_method(
    coefficients,
    *,
    ham_current,
    ham_previous,
    hvib_current,
    hvib_previous,
    dt,
    method=0,
    projector=None,
    rep: Representation = "adiabatic",
    overlap_current=None,
    backend=backend_default,
    propagator: Propagator = exp_propagator,
    nsubsteps: int = 1,
):
    """Propagate amplitudes using the ``Dynamics.cpp`` method selectors.

    Parameters
    ----------
    coefficients
        Batch of electronic coefficient vectors in the old/raw basis.
    ham_current, ham_previous
        Electronic Hamiltonians at ``t+dt`` and ``t``.
    hvib_current, hvib_previous
        Vibronic Hamiltonians at ``t+dt`` and ``t``.
    projector
        ``T_new`` from the C++ code. Its columns express dynamically
        consistent old-state labels in the new raw basis. Consequently
        ``T_new @ C`` maps old coefficients to new raw coefficients, while
        ``T_new† @ H_new @ T_new`` expresses a new Hamiltonian in the old,
        dynamically consistent labels.
    method
        C++ ``electronic_integrator`` selector. Adiabatic methods 0--8 and
        10--15 are implemented, as are aliases 100--115. Method 9 is absent in
        the C++ dispatcher. Diabatic methods 0--3 and aliases 100--103 are
        implemented using their midpoint/split/non-Hermitian forms.
    nsubsteps
        Number of subdivisions used for continuous time evolution. Two-point
        adiabatic methods linearly interpolate their endpoint Hamiltonians at
        substep midpoints. Discrete old-to-new basis mappings are applied
        exactly once per nuclear step.

    Notes
    -----
    Rotation-based options 10--15 are mathematically equivalent to 0--5 for
    exact matrix exponentials and therefore share those kernels here. The C++
    option 8 reads an uninitialized matrix in its first propagation statement;
    Python implements the documented/intended old-then-new half-step form.
    Option 6 intentionally preserves the implemented C++ behavior: a half-step
    with the old vibronic Hamiltonian and no projector application.

    Known C++ bugs and Python choices
    ---------------------------------
    ``method == 8``
        The first C++ propagation uses local variable ``Hvib`` before it has
        been assigned in that branch. Python cannot reproduce undefined
        memory and instead implements the branch's documented intent: an old
        Hvib half-step followed by a transformed-new-Hvib half-step.
    ``method == 12 or method == 112``
        The C++ condition is written ``method == 12 || method == 12``, so alias
        112 is accidentally unreachable. Python accepts 112 consistently with
        the other 100-series aliases.
    ``method == 6``
        C++ forms ``Hvib_old + Hvib_new`` but does not use that sum; it applies
        only ``exp(-i Hvib_old dt/2)`` and no projector. Python preserves this
        implemented behavior, even though it differs from the nearby comment
        describing a two-point scheme.
    """

    method = int(method)
    nsubsteps = int(nsubsteps)
    if nsubsteps < 1:
        raise ValueError("nsubsteps must be positive")

    if method == -1:
        return coefficients
    base_method = method - 100 if 100 <= method < 200 else method
    if rep == "diabatic":
        return _propagate_diabatic_method(
            coefficients,
            ham_current,
            ham_previous,
            hvib_current,
            hvib_previous,
            dt,
            base_method,
            overlap_current,
            backend,
            propagator,
            nsubsteps,
        )
    if rep != "adiabatic":
        raise ValueError("rep must be 'adiabatic' or 'diabatic'")
    if base_method in (10, 11, 12, 13, 14, 15):
        base_method -= 10
    if base_method == 9 or base_method < 0 or base_method > 8:
        raise ValueError(f"unsupported adiabatic electronic integrator {method}")

    current_h = _as_matrix_batch(ham_current, "ham_current")
    previous_h = _as_matrix_batch(ham_previous, "ham_previous")
    current_v = _as_matrix_batch(hvib_current, "hvib_current")
    previous_v = _as_matrix_batch(hvib_previous, "hvib_previous")
    transform = _projector_batch(projector, current_h)
    transform_h = backend.conjugate_transpose(transform)
    current_h_tracked = _rotate_matrix(
        current_h, transform, transform_h, backend
    )
    current_v_tracked = _rotate_matrix(
        current_v, transform, transform_h, backend
    )

    # Endpoint interpolation must occur in one common basis. Coefficients
    # remain in the old dynamically consistent basis during the continuous
    # propagation and are mapped to the current raw basis only at the end.
    if nsubsteps > 1 and base_method in (0, 1, 2):
        state = interpolated_hamiltonian_propagator(
            coefficients,
            previous_h,
            current_h_tracked,
            dt,
            nsubsteps,
            backend,
            propagator,
        )
        return _apply_matrix_to_state(transform, state, backend)
    if nsubsteps > 1 and base_method in (4, 5):
        state = interpolated_hamiltonian_propagator(
            coefficients,
            previous_v,
            current_v_tracked,
            dt,
            nsubsteps,
            backend,
            propagator,
        )
        return _apply_matrix_to_state(transform, state, backend)

    if base_method == 0:
        state = _propagate_fixed_hamiltonian(
            coefficients, previous_h, 0.5 * dt, nsubsteps, backend, propagator
        )
        state = _apply_matrix_to_state(transform, state, backend)
        return _propagate_fixed_hamiltonian(
            state, current_h, 0.5 * dt, nsubsteps, backend, propagator
        )
    if base_method == 1:
        state = _propagate_fixed_hamiltonian(
            coefficients, previous_h, 0.25 * dt, nsubsteps, backend, propagator
        )
        state = _apply_matrix_to_state(transform, state, backend)
        state = _propagate_fixed_hamiltonian(
            state, current_h, 0.5 * dt, nsubsteps, backend, propagator
        )
        state = _apply_matrix_to_state(transform_h, state, backend)
        state = _propagate_fixed_hamiltonian(
            state, previous_h, 0.25 * dt, nsubsteps, backend, propagator
        )
        return _apply_matrix_to_state(transform, state, backend)
    if base_method == 2:
        effective = previous_h + current_h_tracked
        state = _propagate_fixed_hamiltonian(
            coefficients, effective, 0.5 * dt, nsubsteps, backend, propagator
        )
        return _apply_matrix_to_state(transform, state, backend)
    if base_method == 3:
        state = _propagate_fixed_hamiltonian(
            coefficients, previous_v, dt, nsubsteps, backend, propagator
        )
        return _apply_matrix_to_state(transform, state, backend)
    if base_method == 4:
        state = _propagate_fixed_hamiltonian(
            coefficients, previous_v + current_v, 0.5 * dt,
            nsubsteps, backend, propagator,
        )
        return _apply_matrix_to_state(transform, state, backend)
    if base_method == 5:
        effective = previous_v + current_v_tracked
        state = _propagate_fixed_hamiltonian(
            coefficients, effective, 0.5 * dt, nsubsteps, backend, propagator
        )
        return _apply_matrix_to_state(transform, state, backend)
    if base_method == 6:
        return _propagate_fixed_hamiltonian(
            coefficients, previous_v, 0.5 * dt, nsubsteps, backend, propagator
        )
    if base_method == 7:
        state = _propagate_fixed_hamiltonian(
            coefficients, previous_v, 0.5 * dt, nsubsteps, backend, propagator
        )
        state = _apply_matrix_to_state(transform_h, state, backend)
        return _propagate_fixed_hamiltonian(
            state, current_v_tracked, 0.5 * dt,
            nsubsteps, backend, propagator,
        )
    # Intended form of the C++ "new LD" option 8.
    state = _propagate_fixed_hamiltonian(
        coefficients, previous_v, 0.5 * dt, nsubsteps, backend, propagator
    )
    return _propagate_fixed_hamiltonian(
        state, current_v_tracked, 0.5 * dt,
        nsubsteps, backend, propagator,
    )


def _propagate_diabatic_method(
    coefficients,
    ham_current,
    ham_previous,
    hvib_current,
    hvib_previous,
    dt,
    method,
    overlap_current,
    backend,
    propagator,
    nsubsteps,
):
    """Diabatic amplitude branches corresponding to C++ options 0--3."""

    del ham_current, ham_previous
    current = _as_matrix_batch(hvib_current, "hvib_current")
    previous = _as_matrix_batch(hvib_previous, "hvib_previous")
    midpoint = 0.5 * (current + previous)
    if method in (0, 1):
        return _propagate_fixed_hamiltonian(
            coefficients, midpoint, dt, nsubsteps, backend, propagator
        )
    if method == 2:
        return split_step_propagator(
            coefficients, previous, current, dt, backend
        )
    if method == 3:
        if overlap_current is None:
            raise ValueError("diabatic method 3 requires overlap_current")
        overlap = _as_matrix_batch(overlap_current, "overlap_current")
        effective = backend.solve(overlap, midpoint)
        return _propagate_fixed_hamiltonian(
            coefficients, effective, dt, nsubsteps, backend, propagator
        )
    raise ValueError(f"unsupported diabatic electronic integrator {method}")


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
    method: Optional[int] = None,
    nsubsteps: int = 1,
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

    method :
        Optional C++ ``electronic_integrator`` selector. When supplied, the
        old/current Hamiltonians and ``T`` are dispatched through
        :func:`propagate_electronic_method`. When omitted, the original
        one-matrix propagation behavior is retained.

    nsubsteps :
        Number of subdivisions for continuous electronic evolution. Basis
        transformations remain single operations over the full nuclear step.
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

    if method is not None:
        if previous_state is None:
            raise ValueError("method-based propagation requires previous_state")
        C_new = propagate_electronic_method(
            C,
            ham_current=_storage_matrix(storage, traj, rep, "hamiltonian"),
            ham_previous=_snapshot_matrix(previous_state, rep, "hamiltonian"),
            hvib_current=_storage_matrix(storage, traj, rep, "vibronic"),
            hvib_previous=_snapshot_matrix(previous_state, rep, "vibronic"),
            dt=dt,
            method=method,
            projector=T,
            rep=rep,
            overlap_current=(
                storage.ovlp_dia[traj.id, idx] if rep == "diabatic" else None
            ),
            backend=backend,
            propagator=propagator,
            nsubsteps=nsubsteps,
        )
    else:
        # This is a generic basis rotation, not the C++ T_new propagation
        # convention. Method-based dynamics should use the branch above.
        if T is not None:
            C, H = apply_local_diabatization(C, H, T, backend)
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


def _storage_matrix(storage, traj, rep, kind):
    suffix = "adi" if rep == "adiabatic" else "dia"
    prefix = "ham" if kind == "hamiltonian" else "hvib"
    return getattr(storage, f"{prefix}_{suffix}")[traj.id, traj.tbf_ids]


def _snapshot_matrix(state, rep, kind):
    suffix = "adi" if rep == "adiabatic" else "dia"
    prefix = "ham" if kind == "hamiltonian" else "hvib"
    name = f"{prefix}_{suffix}"
    if isinstance(state, dict):
        value = state.get(name)
    else:
        value = getattr(state, name, None)
    if value is None:
        raise ValueError(f"previous_state does not contain {name}")
    return value


def _as_matrix_batch(value, name):
    matrix = value
    if getattr(matrix, "ndim", 0) == 2:
        matrix = matrix[None, ...]
    if getattr(matrix, "ndim", 0) != 3 or matrix.shape[-2] != matrix.shape[-1]:
        raise ValueError(f"{name} must be a square matrix or matrix batch")
    return matrix


def _projector_batch(projector, reference):
    if projector is None:
        nstates = reference.shape[-1]
        identity = backend_default.eye(nstates, dtype=complex)
        return identity[None, ...] if reference.shape[0] == 1 else np.broadcast_to(
            identity, reference.shape
        ).copy()
    transform = _as_matrix_batch(projector, "projector")
    if transform.shape[0] == 1 and reference.shape[0] != 1:
        transform = np.broadcast_to(transform, reference.shape).copy()
    if transform.shape != reference.shape:
        raise ValueError("projector and Hamiltonian batches must have equal shapes")
    return transform
