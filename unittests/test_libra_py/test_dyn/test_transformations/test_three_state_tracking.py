"""Three-state tests for state identity and projector transformations."""

import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.transformations.projectors import (
    compute_permutations,
    compute_projector,
    update_projection,
    update_projectors,
)
from libra_py.dyn.transformations.state_tracking import (
    compute_force_cost_matrix,
    compute_force_cost_matrix_dof_resolved,
    get_reordering,
    get_stochastic_reordering,
    get_stochastic_reordering2,
    get_stochastic_reordering3,
    hungarian_algorithm,
    permutation_matrix,
    permute_states,
)


PERMUTATION = np.array([1, 2, 0])
ENERGIES = np.diag([-0.2, 0.0, 0.3])
CYCLIC_OVERLAP = np.array(
    [
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
    ],
    dtype=complex,
)


def test_three_state_deterministic_algorithms_find_cyclic_permutation():
    np.testing.assert_array_equal(get_reordering(CYCLIC_OVERLAP), PERMUTATION)
    np.testing.assert_array_equal(
        hungarian_algorithm(CYCLIC_OVERLAP, ENERGIES), PERMUTATION
    )


def test_three_state_stochastic_algorithms_are_exact_for_unit_probabilities():
    for function in (
        get_stochastic_reordering,
        get_stochastic_reordering2,
        get_stochastic_reordering3,
    ):
        np.testing.assert_array_equal(
            function(CYCLIC_OVERLAP, np.random.default_rng(9)), PERMUTATION
        )


def test_three_state_projector_corrects_permutation_and_complex_phases():
    phases = np.array([1.0j, -1.0, np.exp(0.25j)])
    overlap = np.zeros((3, 3), dtype=complex)
    overlap[np.arange(3), PERMUTATION] = phases
    params = DynControlParams(
        state_tracking_algo=21,
        do_phase_correction=1,
        phase_correction_tol=1.0e-12,
    )

    projector = compute_projector(params, ENERGIES, overlap)

    np.testing.assert_allclose(overlap @ projector, np.eye(3), atol=1.0e-14)
    np.testing.assert_allclose(projector.conj().T @ projector, np.eye(3))


def test_three_state_batched_permutations_projectors_and_active_states():
    overlaps = np.asarray([CYCLIC_OVERLAP, np.eye(3)])
    energies = np.asarray([ENERGIES, ENERGIES])
    params = DynControlParams(state_tracking_algo=2)

    permutations = compute_permutations(params, energies, overlaps)
    np.testing.assert_array_equal(permutations, [PERMUTATION, [0, 1, 2]])
    np.testing.assert_array_equal(permute_states(permutations, [2, 1]), [0, 1])

    initial = np.asarray([np.eye(3), permutation_matrix([2, 0, 1])])
    updated = update_projectors(params, initial, energies, overlaps)
    np.testing.assert_allclose(updated[0], permutation_matrix(PERMUTATION))
    np.testing.assert_allclose(updated[1], initial[1])


def test_three_state_force_costs_have_state_and_dof_dimensions():
    forces_previous = np.array(
        [[0.0, 0.2], [0.8, -0.1], [-0.4, 0.5]], dtype=float
    )
    forces_current = np.array(
        [[0.1, 0.3], [1.0, -0.2], [-0.2, 0.7]], dtype=float
    )
    arguments = (
        forces_current,
        forces_previous,
        ENERGIES,
        ENERGIES,
        [0.5, -0.2],
        [1.0, 0.5],
        0.25,
    )

    aggregate = compute_force_cost_matrix(*arguments)
    resolved = compute_force_cost_matrix_dof_resolved(*arguments)

    assert aggregate.shape == (3, 3)
    assert resolved.shape == (2, 3, 3)
    np.testing.assert_allclose(np.diag(aggregate), 0.5)
    np.testing.assert_allclose(np.diagonal(resolved, axis1=1, axis2=2), 0.0)


def test_three_state_ld_and_svd_updates_remain_unitary():
    overlap = np.array(
        [
            [0.91, 0.15, -0.02],
            [-0.11, 0.86, 0.18],
            [0.04, -0.14, 0.88],
        ],
        dtype=complex,
    )
    for algorithm in (-1, 5, 6):
        projector = update_projection(
            {"state_tracking_algo": algorithm}, overlap, energies_current=ENERGIES
        )
        np.testing.assert_allclose(
            projector.conj().T @ projector, np.eye(3), atol=1.0e-12
        )
