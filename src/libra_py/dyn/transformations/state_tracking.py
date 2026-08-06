"""State tracking, assignment, phase, and force-cost algorithms.

This module translates the numerical routines in ``src/dyn/dyn_projectors.cpp``.
For a time overlap ``S[i, j] = <psi_i(t)|psi_j(t+dt)>``, row ``i`` is an old
state and column ``j`` is a raw new state. A permutation follows the convention
``perm[i] = j``: old state ``i`` is labelled ``j`` at the new time.

``permutation_matrix(perm)`` returns ``P`` with ``P[perm[i], i] = 1``. Hence
``S @ P`` reorders new-state columns into old-state order, and column ``i`` of
``P`` identifies the raw new label of active old state ``i``.

Known C++ differences
---------------------
``dyn_projectors.cpp::get_reordering`` composes cumulative transpositions in
an order that returns the inverse mapping for cycles longer than two. The
Python function keeps the documented old-to-new convention required by
``permutation2cmatrix``, active-state reassignment, and electronic propagation.

The C++ DOF-resolved force-cost routine uses the ``i``-state acceleration when
constructing both the ``i`` and ``j`` trial displacements. Python presently
preserves that expression for numerical compatibility and calls it out in the
function docstring rather than silently changing the algorithm.
"""

from __future__ import annotations

from itertools import permutations

import numpy as np


def compute_phase_corrections(time_overlap, tol=1.0e-3):
    """Return phases ``S_ii / |S_ii|``, or one below the overlap threshold.

    ``time_overlap`` should already be permuted so that matched states are on
    its diagonal. For diagonal element ``z_i`` above ``tol``, the returned
    factor is ``f_i = z_i / |z_i|``; otherwise it is one. A new eigenvector is
    corrected by ``conj(f_i)``, so the corrected diagonal overlap becomes
    ``|z_i|``. The return value has shape ``(nstates,)`` and unit modulus.
    """

    overlap = _square(time_overlap, "time_overlap", complex)
    diagonal = np.diag(overlap)
    phases = np.ones(len(diagonal), dtype=complex)
    mask = np.abs(diagonal) > float(tol)
    phases[mask] = diagonal[mask] / np.abs(diagonal[mask])
    return phases


def get_reordering(time_overlap):
    """Greedily reorder columns by maximum overlap (algorithm 1).

    The result follows the public C++ convention ``perm[old] = new`` so that
    ``time_overlap @ permutation_matrix(perm)`` moves selected overlaps onto
    the diagonal. The C++ cumulative-swap implementation composes in the
    opposite order and therefore returns the inverse for cycles longer than
    two; Python corrects that inconsistency because projectors, active-state
    reassignment, and ``Dynamics.cpp`` propagation all require old-to-new.

    .. warning::
       This intentionally differs from the current C++ implementation for
       permutations containing cycles of length three or greater. The C++
       result is the inverse of its own documented ``perm[old] = new`` mapping
       in those cases. Two-state tests cannot reveal this bug.
    """

    work = np.array(_square(time_overlap, "time_overlap", complex), copy=True)
    nstates = len(work)
    cumulative = np.arange(nstates)
    for column in range(nstates):
        for _ in range(nstates):
            target = int(np.argmax(np.abs(work[:, column])))
            if target == column:
                break
            work[:, [column, target]] = work[:, [target, column]]
            swap = np.arange(nstates)
            swap[[column, target]] = swap[[target, column]]
            # Compose on the right to retain perm[old] = new. This differs
            # from the buggy C++ cycle composition but matches all consumers.
            cumulative = cumulative[swap]
    return cumulative


def make_cost_matrix(overlaps, energies, alpha=0.0, scaling_function=0):
    """Build the energy-aware overlap assignment score ``|S_ab|² f(dE)``.

    The selectors reproduce the *implemented* C++ branches:

    ``0`` or other
        no energy scaling;
    ``1``
        ``exp[-(alpha*dE)^2]``;
    ``2``
        ``exp[-alpha*|dE|]``;
    ``3``
        ``exp[-alpha*max(dE, 0)]``.
    """

    overlap = _square(overlaps, "overlaps", complex)
    energy = _energy_diagonal(energies, len(overlap))
    gaps = energy[:, None] - energy[None, :]
    if scaling_function == 1:
        scale = np.exp(-(float(alpha) * gaps) ** 2)
    elif scaling_function == 2:
        scale = np.exp(-float(alpha) * np.abs(gaps))
    elif scaling_function == 3:
        scale = np.exp(-float(alpha) * np.maximum(gaps, 0.0))
    else:
        scale = 1.0
    return np.abs(overlap) ** 2 * scale


# C++-style spelling retained as an alias.
make_cost_mat = make_cost_matrix


def hungarian_algorithm(cost_matrix, energies=None, alpha=0.0, scaling_function=0):
    """Return the minimum-cost assignment, or maximize overlap scores.

    With ``energies`` supplied, ``cost_matrix`` is interpreted as a time
    overlap and the energy-aware score is maximized, matching the C++ overload.
    """

    if energies is not None:
        costs = -make_cost_matrix(cost_matrix, energies, alpha, scaling_function)
    else:
        costs = np.asarray(cost_matrix, dtype=float)
    costs = _square(costs, "cost_matrix", float)
    return _linear_assignment(costs)


def munkres_kuhn(overlaps, energies, alpha=0.0, verbosity=0, scaling_function=0):
    """Maximize the energy-aware overlap score (C++ algorithm 2)."""

    del verbosity
    return hungarian_algorithm(overlaps, energies, alpha, scaling_function)


# C++ public name.
Munkres_Kuhn = munkres_kuhn


def get_stochastic_reordering(time_overlap, rng=None):
    """Sequential stochastic assignment without replacement (algorithm 3)."""

    probabilities = np.abs(_square(time_overlap, "time_overlap", complex)) ** 2
    rng = _rng(rng)
    available = list(range(len(probabilities)))
    result = np.empty(len(probabilities), dtype=int)
    for source in range(len(probabilities)):
        weights = probabilities[source, available]
        selected = _sample(weights, rng)
        result[source] = available.pop(selected)
    return result


def get_stochastic_reordering2(time_overlap, rng=None):
    """Sample complete permutations by products of overlap probabilities (32)."""

    probabilities = np.abs(_square(time_overlap, "time_overlap", complex)) ** 2
    choices = list(permutations(range(len(probabilities))))
    weights = np.asarray([
        np.prod(probabilities[np.arange(len(choice)), choice]) for choice in choices
    ])
    return np.asarray(choices[_sample(weights, _rng(rng))], dtype=int)


def get_stochastic_reordering3(
    time_overlap,
    rng=None,
    convergence=0,
    max_number_of_attempts=100,
    filter_tol=0.0,
    verbosity_level=0,
):
    """Independently sample rows until a valid permutation is found (33).

    If no permutation is found, ``convergence=0`` returns the identity while
    ``convergence=1`` raises ``RuntimeError`` (the safe Python analogue of the
    C++ process exit).
    """

    probabilities = np.abs(_square(time_overlap, "time_overlap", complex)) ** 2
    probabilities[probabilities <= float(filter_tol)] = 0.0
    rng = _rng(rng)
    for attempt in range(int(max_number_of_attempts)):
        chosen = np.asarray([_sample(row, rng) for row in probabilities], dtype=int)
        if verbosity_level:
            print(f"state-reordering attempt {attempt}: {chosen.tolist()}")
        if len(np.unique(chosen)) == len(chosen):
            return chosen
    if convergence:
        raise RuntimeError(
            f"stochastic reordering did not converge in {max_number_of_attempts} attempts"
        )
    return np.arange(len(probabilities))


def permutation_matrix(permutation):
    """Convert an old-to-new state mapping into a projector matrix.

    For ``perm[i] = j``, set ``P[j, i] = 1``. Therefore column ``i`` is old
    tracked state ``i`` expressed in the raw new-state basis, ``S @ P`` puts
    matched overlaps on the diagonal, and ``P @ C_old`` gives raw new-basis
    coefficients. The result is a complex unitary matrix.
    """

    permutation = _permutation(permutation)
    result = np.zeros((len(permutation), len(permutation)), dtype=complex)
    result[permutation, np.arange(len(permutation))] = 1.0
    return result


permutation2cmatrix = permutation_matrix


def permute_states(permutations_by_trajectory, active_states):
    """Map every active state through its trajectory-specific permutation.

    ``permutations_by_trajectory[t, i]`` is the new label of old state ``i``
    on trajectory ``t``. The result therefore contains
    ``permutations_by_trajectory[t, active_states[t]]``.
    """

    perms = np.asarray(permutations_by_trajectory, dtype=int)
    states = np.asarray(active_states, dtype=int)
    if perms.ndim != 2 or states.shape != (perms.shape[0],):
        raise ValueError("permutations and active_states dimensions do not agree")
    return perms[np.arange(len(states)), states]


def compute_force_cost_matrix(
    forces_current, forces_previous, energies_current, energies_previous,
    momentum, inverse_mass, dt, active_state=0,
):
    """Detect predicted pairwise energy-gap sign changes from forces.

    This is ``compute_F_cost_matrix``. ``energies_current`` and
    ``active_state`` are accepted for signature compatibility; the C++ formula
    presently uses the previous energies and does not use the active state.
    """

    del energies_current, active_state
    fc, fp, ep, p, im = _force_inputs(
        forces_current, forces_previous, energies_previous, momentum, inverse_mass
    )
    acceleration = fc * im[None, :]
    displacement = p[None, :] * im[None, :] * dt + 0.5 * acceleration * dt**2
    energy_change = -np.sum(fp * displacement, axis=1)
    gap = ep[None, :] - ep[:, None]
    new_gap = gap + energy_change[None, :] - energy_change[:, None]
    return 0.5 * (1.0 - np.sign(gap) * np.sign(new_gap)) + 0.0j


compute_F_cost_matrix = compute_force_cost_matrix


def compute_force_cost_matrix_dof_resolved(
    forces_current, forces_previous, energies_current, energies_previous,
    momentum, inverse_mass, dt, active_state=0,
):
    """Return the C++ force-gap sign product separately for every DOF.

    The result has shape ``(ndof, nstates, nstates)``. Values ``-1`` indicate
    a predicted sign change of the corresponding pairwise energy gap, ``1``
    indicates no sign change, and diagonal/zero-gap entries are zero.

    .. warning::
       The C++ source computes ``dq_j`` with acceleration ``a_i`` rather than
       ``a_j``. This is likely a typo because the aggregate
       ``compute_F_cost_matrix`` uses the separate ``j``-state acceleration.
       Python preserves the C++ expression here so translated calculations
       remain comparable; the behavior is covered by regression tests.
    """

    del energies_current, active_state
    fc, fp, ep, p, im = _force_inputs(
        forces_current, forces_previous, energies_previous, momentum, inverse_mass
    )
    result = np.empty((fc.shape[1], fc.shape[0], fc.shape[0]))
    gap = ep[None, :] - ep[:, None]
    for dof in range(fc.shape[1]):
        # Preserve the C++ expression, including its use of the i-state
        # acceleration for both i and j displacements.
        for i in range(fc.shape[0]):
            dq_i = p[dof] * im[dof] * dt + 0.5 * fc[i, dof] * im[dof] * dt**2
            de_i = -fp[i, dof] * dq_i
            de_j = -fp[:, dof] * dq_i
            new_gap = gap[i] + de_j - de_i
            result[dof, i] = np.sign(gap[i]) * np.sign(new_gap)
    return result


compute_F_cost_matrix_dof_resolved = compute_force_cost_matrix_dof_resolved


def compute_force_cost_matrix2(forces_current, forces_previous, *unused):
    """Experimental normalized force-overlap score from the C++ placeholder."""

    del unused
    current = _forces(forces_current, "forces_current")
    previous = _forces(forces_previous, "forces_previous")
    if current.shape != previous.shape:
        raise ValueError("current and previous force arrays must have equal shape")
    result = np.empty((len(current), len(current)), dtype=complex)
    for i in range(len(current)):
        for j in range(len(current)):
            ni = np.dot(previous[i], previous[i])
            nj = np.dot(current[j], current[j])
            result[i, j] = (
                np.dot(previous[i], current[j]) / np.sqrt(ni * nj)
                if ni >= 1.0e-12 and nj >= 1.0e-12 else float(i == j)
            )
    return result


compute_F_cost_matrix2 = compute_force_cost_matrix2


def _linear_assignment(costs):
    """O(N^3) square Hungarian minimization returning row-to-column labels."""

    n = len(costs)
    u = np.zeros(n + 1)
    v = np.zeros(n + 1)
    match = np.zeros(n + 1, dtype=int)
    path = np.zeros(n + 1, dtype=int)
    for row in range(1, n + 1):
        match[0] = row
        col0 = 0
        minimum = np.full(n + 1, np.inf)
        used = np.zeros(n + 1, dtype=bool)
        while True:
            used[col0] = True
            row0 = match[col0]
            delta = np.inf
            col1 = 0
            for col in range(1, n + 1):
                if not used[col]:
                    current = costs[row0 - 1, col - 1] - u[row0] - v[col]
                    if current < minimum[col]:
                        minimum[col], path[col] = current, col0
                    if minimum[col] < delta:
                        delta, col1 = minimum[col], col
            for col in range(n + 1):
                if used[col]:
                    u[match[col]] += delta
                    v[col] -= delta
                else:
                    minimum[col] -= delta
            col0 = col1
            if match[col0] == 0:
                break
        while True:
            col1 = path[col0]
            match[col0] = match[col1]
            col0 = col1
            if col0 == 0:
                break
    assignment = np.empty(n, dtype=int)
    for col in range(1, n + 1):
        assignment[match[col] - 1] = col - 1
    return assignment


def _sample(weights, rng):
    weights = np.asarray(weights, dtype=float)
    total = weights.sum()
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("state-tracking probabilities have zero or invalid norm")
    value = float(rng.uniform(0.0, 1.0)) * total
    return min(int(np.searchsorted(np.cumsum(weights), value, side="right")), len(weights) - 1)


def _rng(rng):
    return np.random.default_rng() if rng is None else rng


def _square(value, name, dtype):
    array = np.asarray(value, dtype=dtype)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    return array


def _energy_diagonal(value, nstates):
    array = np.asarray(value)
    diagonal = np.real(np.diag(array)) if array.ndim == 2 else np.real(array).reshape(-1)
    if diagonal.shape != (nstates,):
        raise ValueError("energies must be a state vector or square matrix")
    return diagonal


def _permutation(value):
    perm = np.asarray(value, dtype=int)
    if perm.ndim != 1 or sorted(perm.tolist()) != list(range(len(perm))):
        raise ValueError("permutation must contain every index exactly once")
    return perm


def _forces(value, name):
    array = np.real(np.asarray(value))
    if array.ndim != 2:
        raise ValueError(f"{name} must have shape (nstates, ndof)")
    return array


def _force_inputs(fc, fp, ep, momentum, inverse_mass):
    fc = _forces(fc, "forces_current")
    fp = _forces(fp, "forces_previous")
    if fc.shape != fp.shape:
        raise ValueError("current and previous force arrays must have equal shape")
    ep = _energy_diagonal(ep, len(fc))
    p = np.asarray(momentum, dtype=float).reshape(-1)
    im = np.asarray(inverse_mass, dtype=float).reshape(-1)
    if p.shape != im.shape or len(p) != fc.shape[1]:
        raise ValueError("momentum and inverse_mass must match the force DOFs")
    return fc, fp, ep, p, im
