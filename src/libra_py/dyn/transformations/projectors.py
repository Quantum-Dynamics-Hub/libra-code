"""Projection construction and update workflows translated from C++.

Conventions
-----------
Let ``S[i, j] = <psi_i(t) | psi_j(t+dt)>``. Rows therefore carry old state
labels and columns carry the raw new labels. A state permutation is stored as
``perm[i] = j``: old state ``i`` has raw new label ``j``.

The corresponding projector ``P`` is defined by ``P[perm[i], i] = 1``. Its
columns are the new raw states arranged in the old, dynamically consistent
order. Consequently:

* ``S @ P`` puts matched overlaps on the diagonal;
* ``argmax(abs(P[:, i]))`` reassigns active old state ``i``;
* ``P @ C_old`` maps coefficients into the new raw basis;
* ``P.conj().T @ H_new @ P`` expresses a new raw Hamiltonian using the old
  dynamically consistent labels.

Phase correction right-multiplies ``P`` by ``diag(conj(f))``, where
``f_i = (S @ P)_ii / |(S @ P)_ii|``. Thus the diagonal of the corrected
``S @ P`` is real and nonnegative. These are the conventions used by
``dyn_projectors.cpp``, ``dyn_ham.cpp``, and the projector-aware branches of
``Dynamics.cpp::propagate_electronic``.

Compatibility note
------------------
The Python greedy tracker deliberately fixes the multi-state cycle-composition
bug in ``dyn_projectors.cpp::get_reordering``. Thus algorithm 1 can differ from
the current C++ result for cycles of length three or more, while algorithms 2
and 21 already follow the projector convention directly. See
``state_tracking.get_reordering`` for the precise difference.
"""

from __future__ import annotations

import numpy as np

from .local_diabatization import orthogonalized_T
from .state_tracking import (
    compute_force_cost_matrix,
    compute_phase_corrections,
    get_reordering,
    get_stochastic_reordering,
    get_stochastic_reordering2,
    get_stochastic_reordering3,
    hungarian_algorithm,
    munkres_kuhn,
    permutation_matrix,
)


def compute_projector(params, energies, time_overlap, rng=None):
    """Compute one instantaneous state-tracking and phase projector.

    Parameters
    ----------
    params
        ``DynControlParams`` or dictionary. Relevant keys are
        ``state_tracking_algo``, ``MK_alpha``, ``MK_scaling_function``,
        ``do_phase_correction``, and ``phase_correction_tol``.
    energies
        Current state energies, either ``(nstates,)`` or a square Hamiltonian.
        They enter energy-aware assignment algorithms 2, 21, and 4.
    time_overlap
        Square matrix ``S[old, new]``.
    rng
        Random-number generator used by stochastic algorithms 3, 32, and 33.

    Returns
    -------
    numpy.ndarray
        Unitary permutation/phase matrix ``P`` with shape
        ``(nstates, nstates)``. Before phase scaling,
        ``P[perm[i], i] = 1``; after scaling the location is unchanged and only
        its unit-modulus value changes.

    Notes
    -----
    State-tracking option numbers are the C++ numbers: 1 greedy, 2/4
    Munkres-Kuhn, 21 Hungarian, and 3/32/33 stochastic variants. If phase
    correction is enabled, phases are computed only after permutation from
    ``S @ P`` and are applied to projector columns.
    """

    overlap = _square(time_overlap, "time_overlap")
    algorithm = int(_param(params, "state_tracking_algo", 0))
    permutation = _tracking_permutation(algorithm, params, overlap, energies, rng)
    projector = permutation_matrix(permutation)
    reordered_overlap = overlap @ projector
    if int(_param(params, "do_phase_correction", 0)):
        phases = compute_phase_corrections(
            reordered_overlap, float(_param(params, "phase_correction_tol", 1.0e-3))
        )
        projector = projector @ np.diag(np.conjugate(phases))
    return projector


def compute_permutations(params, energies, time_overlaps, rng=None):
    """Compute state permutations for a batch of trajectories.

    Parameters are read from either ``DynControlParams`` or a dictionary.
    ``energies`` can be a shared state vector/matrix or a trajectory batch;
    ``time_overlaps`` has shape ``(ntraj, nstates, nstates)``. The returned
    integer array has one old-to-new mapping per processed trajectory. With
    ``isNBRA=1``, only the first trajectory is processed, matching C++.
    """

    overlaps = _batch_square(time_overlaps, "time_overlaps")
    energy_batch = _batch_energies(energies, len(overlaps), overlaps.shape[-1])
    count = 1 if int(_param(params, "isNBRA", 0)) else len(overlaps)
    return np.asarray([
        _tracking_permutation(
            int(_param(params, "state_tracking_algo", 0)),
            params, overlaps[index], energy_batch[index], rng,
        )
        for index in range(count)
    ])


def update_projectors(params, projectors, energies, time_overlaps, rng=None):
    """Update cumulative projectors as C++ ``update_projectors`` does.

    ``projectors[t]`` is the already accumulated raw-to-dynamically-consistent
    basis map for trajectory ``t``. For instantaneous permutation ``p_t``, the
    update is

    ``P_new = P_old @ p_t``.

    Phase factors are evaluated from

    ``S_consistent = P_old† @ S_raw @ P_new``

    and then applied to columns of ``P_new``. ``projectors`` and
    ``time_overlaps`` both have shape ``(ntraj, nstates, nstates)``; the return
    value has the same shape and the inputs are not modified.
    """

    old = _batch_square(projectors, "projectors")
    overlaps = _batch_square(time_overlaps, "time_overlaps")
    if old.shape != overlaps.shape:
        raise ValueError("projectors and time_overlaps must have equal shapes")
    energy_batch = _batch_energies(energies, len(old), old.shape[-1])
    result = np.array(old, copy=True)
    count = 1 if int(_param(params, "isNBRA", 0)) else len(old)
    for index in range(count):
        instantaneous = compute_projector(
            {**_params_dict(params), "do_phase_correction": 0},
            energy_batch[index], overlaps[index], rng,
        )
        result[index] = old[index] @ instantaneous
        transformed_overlap = old[index].conj().T @ overlaps[index] @ result[index]
        if int(_param(params, "do_phase_correction", 0)):
            phases = compute_phase_corrections(
                transformed_overlap,
                float(_param(params, "phase_correction_tol", 1.0e-3)),
            )
            result[index] = result[index] @ np.diag(np.conjugate(phases))
    return result


def update_projection(
    params,
    time_overlap,
    *,
    current_projector=None,
    energies_current=None,
    energies_previous=None,
    forces_current=None,
    forces_previous=None,
    momentum=None,
    inverse_mass=None,
    dt=None,
    active_state=0,
    rng=None,
):
    """Compute one ``update_proj_adi`` transformation from array inputs.

    This is the array-level equivalent of one trajectory iteration in
    ``dyn_ham.cpp::update_proj_adi``. ``time_overlap`` always uses
    ``S[old, new]`` orientation, and the returned ``T_new`` uses the projector
    convention described in this module's docstring.

    Algorithms
    ----------
    ``-1``
        Local diabatization: polar-orthogonalize ``inv(S)`` using
        :func:`orthogonalized_T` from ``local_diabatization.py``.
    ``0``
        Keep ``current_projector`` unchanged.
    ``1, 2, 21, 3, 32, 33``
        Build an instantaneous permutation/phase projector.
    ``4``
        Construct a force-based crossing cost matrix, then use the option-2
        assignment and phase workflow. Current/previous energies and forces,
        momentum, inverse mass, and ``dt`` are required.
    ``5``
        SVD polar factor ``U @ V†`` of the overlap.
    ``6``
        Adaptive regularized SVD local-diabatization transformation.

    Returns
    -------
    numpy.ndarray
        Square ``T_new`` used by electronic propagation as ``T_new @ C`` and
        by corrected Hamiltonian branches as ``T_new† @ H_new @ T_new``.
    """

    overlap = _square(time_overlap, "time_overlap")
    nstates = len(overlap)
    current = (
        np.eye(nstates, dtype=complex)
        if current_projector is None else _square(current_projector, "current_projector")
    )
    algorithm = int(_param(params, "state_tracking_algo", -1))
    if algorithm == -1:
        try:
            inverse = np.linalg.inv(overlap)
        except np.linalg.LinAlgError as exc:
            raise ValueError("time-overlap matrix is singular") from exc
        return orthogonalized_T(inverse)
    if algorithm == 0:
        return np.array(current, copy=True)
    if algorithm in (1, 2, 21, 3, 32, 33):
        return compute_projector(params, energies_current, overlap, rng)
    if algorithm == 4:
        required = {
            "energies_current": energies_current,
            "energies_previous": energies_previous,
            "forces_current": forces_current,
            "forces_previous": forces_previous,
            "momentum": momentum,
            "inverse_mass": inverse_mass,
            "dt": dt,
        }
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise ValueError("force-based tracking requires " + ", ".join(missing))
        force_cost = compute_force_cost_matrix(
            forces_current, forces_previous, energies_current, energies_previous,
            momentum, inverse_mass, float(dt), active_state,
        )
        return compute_projector(params, energies_current, force_cost, rng)
    if algorithm == 5:
        left, _, right_h = np.linalg.svd(overlap)
        return left @ right_h
    if algorithm == 6:
        left, singular, right_h = np.linalg.svd(overlap)
        right = right_h.conj().T
        regularizer = singular / (singular + 1.0e-12)
        f = np.diag(regularizer)
        identity = np.eye(nstates, dtype=complex)
        projection = right @ f @ right.conj().T
        t_tilde = left @ f @ right_h + identity - projection
        u2, _, vh2 = np.linalg.svd(t_tilde)
        return (u2 @ vh2).conj().T
    raise ValueError(f"unsupported state_tracking_algo {algorithm}")


def update_proj_adi(params, storage, traj, previous=None, rng=None):
    """Update ``storage.proj_adi`` for active TBFs using storage tensors.

    The function reads ``time_overlap_adi``, ``ham_adi``, ``proj_adi``, nuclear
    momenta, and inverse masses for ``traj.tbf_ids``. It calls
    :func:`update_projection` independently for each active TBF and writes the
    resulting batch back to ``storage.proj_adi``.

    ``previous`` may be another storage-like object or a dictionary snapshot
    containing ``ham_adi`` and ``d1ham_adi`` batches. For force-based option 4
    it is required; other algorithms only use current storage data. The
    returned shape is ``(ntbf_active, nstates, nstates)``.

    Downstream use
    --------------
    ``DynamicsEngine`` passes this matrix batch to the method-aware electronic
    propagator. It also maps an active old state ``i`` to the raw new state at
    ``argmax(abs(P[:, i]))``. Phase factors never affect that index mapping.
    """

    idx = traj.tbf_ids
    overlaps = np.asarray(storage.time_overlap_adi[traj.id, idx])
    energies = np.asarray(storage.ham_adi[traj.id, idx])
    current = np.asarray(storage.proj_adi[traj.id, idx])
    momenta = np.asarray(storage.p[traj.id, idx])
    inverse_masses = np.asarray(storage.iM[traj.id, idx])
    algorithm = int(_param(params, "state_tracking_algo", -1))
    previous_energies = _previous_field(previous, "ham_adi", traj, idx)
    previous_derivatives = _previous_field(previous, "d1ham_adi", traj, idx)
    current_derivatives = getattr(storage, "d1ham_adi", None)
    if current_derivatives is not None:
        current_derivatives = np.asarray(current_derivatives[traj.id, idx])

    result = np.empty_like(current, dtype=complex)
    for position in range(len(idx)):
        kwargs = {}
        if algorithm == 4:
            if current_derivatives is None or previous_derivatives is None:
                raise ValueError("force-based tracking requires current and previous d1ham_adi")
            kwargs = {
                "energies_previous": previous_energies[position],
                "forces_current": -np.real(np.diagonal(
                    current_derivatives[position], axis1=-2, axis2=-1
                )).T,
                "forces_previous": -np.real(np.diagonal(
                    previous_derivatives[position], axis1=-2, axis2=-1
                )).T,
                "momentum": momenta[position],
                "inverse_mass": inverse_masses[position],
                "dt": float(_param(params, "dt", 41.0)),
                "active_state": int(storage.act_states[traj.id, idx[position]]),
            }
        result[position] = update_projection(
            params,
            overlaps[position],
            current_projector=current[position],
            energies_current=energies[position],
            rng=rng,
            **kwargs,
        )
    storage.proj_adi[traj.id, idx] = result
    return result


def _tracking_permutation(algorithm, params, overlap, energies, rng):
    if algorithm in (-1, 0, 5, 6):
        return np.arange(len(overlap))
    if algorithm == 1:
        return get_reordering(overlap)
    if algorithm in (2, 4):
        return munkres_kuhn(
            overlap, energies, _param(params, "MK_alpha", 0.0),
            _param(params, "MK_verbosity", 0),
            _param(params, "MK_scaling_function", 0),
        )
    if algorithm == 21:
        return hungarian_algorithm(
            overlap, energies, _param(params, "MK_alpha", 0.0),
            _param(params, "MK_scaling_function", 0),
        )
    if algorithm == 3:
        return get_stochastic_reordering(overlap, rng)
    if algorithm == 32:
        return get_stochastic_reordering2(overlap, rng)
    if algorithm == 33:
        return get_stochastic_reordering3(
            overlap, rng,
            _param(params, "convergence", 0),
            _param(params, "max_number_attempts", 100),
            _param(params, "min_probability_reordering", 0.0),
        )
    raise ValueError(f"unsupported state_tracking_algo {algorithm}")


def _previous_field(previous, name, traj, idx):
    if previous is None:
        return None
    if isinstance(previous, dict):
        value = previous.get(name)
        return None if value is None else np.asarray(value)
    value = getattr(previous, name, None)
    return None if value is None else np.asarray(value[traj.id, idx])


def _batch_energies(value, count, nstates):
    array = np.asarray(value)
    if array.ndim == 1 or (array.ndim == 2 and array.shape == (nstates, nstates)):
        array = np.broadcast_to(array, (count,) + array.shape)
    if len(array) != count:
        raise ValueError("energies batch dimension does not agree")
    return array


def _square(value, name):
    array = np.asarray(value, dtype=complex)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    return array


def _batch_square(value, name):
    array = np.asarray(value, dtype=complex)
    if array.ndim == 2:
        array = array[None, ...]
    if array.ndim != 3 or array.shape[-2] != array.shape[-1]:
        raise ValueError(f"{name} must be a square matrix batch")
    return array


def _param(params, name, default):
    return params.get(name, default) if isinstance(params, dict) else getattr(params, name, default)


def _params_dict(params):
    if isinstance(params, dict):
        return dict(params)
    return dict(vars(params))
