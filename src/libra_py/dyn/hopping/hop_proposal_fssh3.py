"""FSSH3 hop proposals translated from ``dyn_hop_proposal_fssh3.cpp``.

FSSH3 infers a population-transfer matrix from the populations at two
successive time steps, then distributes the observed loss of the active-state
population among positive fluxes leading to other states.  The implementation
below follows the clean, public C++ routine; the experimental ``*_dev`` routine
is intentionally not part of the public API in either implementation.
"""

from __future__ import annotations

import numpy as np


def adjust_signs(matrix):
    """Adjust matrix-column signs using the first row and column.

    For every column ``i > 0``, the C++ helper computes

        s_i = sign(J_i0 J_0i)

    (using ``s_i = 1`` when the product is zero) and scales column ``i`` by
    ``-s_i``.  A copy is returned; the input is not modified.
    """

    result = _square_real(matrix, "matrix").copy()
    source = np.asarray(matrix, dtype=float)
    for i in range(1, result.shape[1]):
        scale = np.sign(source[i, 0] * source[0, i])
        result[:, i] *= -(scale if scale != 0.0 else 1.0)
    return result


def find_best_matrix(c_new, c_old, matrix):
    """Choose the lowest-error sign variant considered by the C++ helper.

    Four candidates are compared: the original matrix, its first column sign
    flipped, its first row sign flipped, and both flipped.  :func:`adjust_signs`
    is applied to each candidate.  The candidate minimizing
    ``||c_new - J c_old||²`` is returned.

    This replaces the C++ ``sum_row`` expression noted as bugged in the source
    with its documented squared-error intent.
    """

    c_new = _column(c_new, "c_new")
    c_old = _column(c_old, "c_old")
    matrix = _square_real(matrix, "matrix")
    if matrix.shape[0] != len(c_new) or len(c_old) != len(c_new):
        raise ValueError("matrix and population-vector dimensions must match")

    variants = []
    for flip_row, flip_column in ((False, False), (False, True), (True, False), (True, True)):
        candidate = matrix.copy()
        if flip_column:
            candidate[:, 0] *= -1.0
        if flip_row:
            candidate[0, :] *= -1.0
        variants.append(adjust_signs(candidate))
    errors = [np.sum((c_new - candidate @ c_old) ** 2) for candidate in variants]
    return variants[int(np.argmin(errors))]


def hopping_probabilities_fssh3(
    params, density, density_old, active_state, errors=None
):
    """Compute FSSH3 proposal probabilities and optimization diagnostics.

    Let ``P_old`` and ``P_new`` be the diagonal populations at successive
    dynamics steps separated by ``dt``.  The C++ routine obtains a matrix ``J``
    by steepest-descent minimization of

        L0 = ||J P_old - y||²,
        dL0/dJ = 2 (J P_old - y) P_oldᵀ.

    For ``fssh3_approach_option == 0``, ``y = P_new`` and ``J`` starts as the
    identity (a transition-matrix interpretation).  For options 1 and 2,
    ``y = (P_new-P_old)/dt`` and ``J`` is a flux/rate matrix initialized to zero
    or to nearest-neighbor ±1 entries, respectively.  Optimization uses
    ``fssh3_dt`` as its gradient step, at most ``fssh3_max_steps`` iterations,
    and stops at ``fssh3_err_tol``.

    The total observed probability of leaving active state ``i`` is

        P_out = clip[-dt (dP_i/dt) / P_old_i, 0, 1].

    Positive entries ``J_ji`` for ``j != i`` determine the relative target
    weights.  They are normalized to sum to ``P_out``; the staying probability
    is ``1-P_out``.  If there is no positive target flux, this Python version
    safely leaves the trajectory on ``i`` instead of reproducing the C++
    division by zero.

    ``errors``, when supplied, is updated in place with the five C++ diagnostic
    values ``[Lagrangian, L0, L1, L2, L3]``.  Only the unconstrained ``L0`` term
    is active in the clean public FSSH3 routine, so the last three values are
    zero.
    """

    density = np.asarray(density, dtype=complex)
    density_old = np.asarray(density_old, dtype=complex)
    if density.ndim != 2 or density.shape[0] != density.shape[1]:
        raise ValueError("density must be a square matrix")
    if density_old.shape != density.shape:
        raise ValueError("density and density_old dimensions must match")
    nstates = density.shape[0]
    active_state = int(active_state)
    if not 0 <= active_state < nstates:
        raise ValueError(f"active_state must be in [0, {nstates})")

    old = np.diag(density_old).real
    new = np.diag(density).real
    dt = float(_param(params, "dt", 41.0))
    if dt == 0.0:
        raise ValueError("dt must be nonzero")
    approach = int(_param(params, "fssh3_approach_option", 0))
    if approach not in (0, 1, 2):
        raise ValueError("fssh3_approach_option must be 0, 1, or 2")

    target = new if approach == 0 else (new - old) / dt
    if approach == 0:
        transfer = np.eye(nstates)
    elif approach == 1:
        transfer = np.zeros((nstates, nstates))
    else:
        transfer = np.zeros((nstates, nstates))
        indices = np.arange(nstates - 1)
        transfer[indices, indices + 1] = 1.0
        transfer[indices + 1, indices] = -1.0

    step = float(_param(params, "fssh3_dt", 0.001))
    max_steps = int(_param(params, "fssh3_max_steps", 1000))
    tolerance = float(_param(params, "fssh3_err_tol", 1.0e-7))
    loss = np.inf
    for _ in range(max_steps):
        residual = transfer @ old - target
        loss = float(residual @ residual)
        if abs(loss) <= tolerance:
            break
        transfer -= step * 2.0 * np.outer(residual, old)
    residual = transfer @ old - target
    loss = float(residual @ residual)
    diagnostics = np.array([loss, loss, 0.0, 0.0, 0.0])
    _store_errors(errors, diagnostics)

    old_active = old[active_state]
    probability_out = 0.0
    if old_active > 0.0:
        probability_out = np.clip(-(new[active_state] - old_active) / old_active, 0.0, 1.0)
    weights = np.maximum(transfer[:, active_state], 0.0)
    weights[active_state] = 0.0
    norm = weights.sum()
    result = np.zeros(nstates, dtype=float)
    if norm > 0.0:
        result = weights * probability_out / norm
        result[active_state] = 1.0 - probability_out
    else:
        result[active_state] = 1.0
    return result


def _store_errors(errors, values):
    if errors is None:
        return
    if len(errors) != 5:
        raise ValueError("errors must contain five elements")
    for index, value in enumerate(values):
        errors[index] = float(value)


def _param(params, name, default):
    return params.get(name, default) if isinstance(params, dict) else getattr(params, name, default)


def _square_real(value, name):
    array = np.asarray(value, dtype=float)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(f"{name} must be a square matrix")
    return array


def _column(value, name):
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    return array
