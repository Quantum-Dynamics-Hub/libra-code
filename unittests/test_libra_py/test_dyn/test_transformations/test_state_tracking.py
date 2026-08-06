import numpy as np
import pytest

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.transformations.projectors import (
    compute_projector,
    update_projection,
)
from libra_py.dyn.transformations.state_tracking import (
    compute_force_cost_matrix,
    compute_phase_corrections,
    get_reordering,
    get_stochastic_reordering2,
    get_stochastic_reordering3,
    hungarian_algorithm,
    make_cost_matrix,
    permutation_matrix,
    permute_states,
)


SWAP = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)


def test_phase_corrections_follow_diagonal_overlap_phases_and_tolerance():
    overlap = np.diag([2.0j, 1.0e-5j, -3.0])
    np.testing.assert_allclose(
        compute_phase_corrections(overlap, tol=1.0e-3),
        [1.0j, 1.0, -1.0],
    )


def test_greedy_and_hungarian_tracking_find_swapped_states():
    energies = np.diag([0.0, 0.1])
    np.testing.assert_array_equal(get_reordering(SWAP), [1, 0])
    np.testing.assert_array_equal(hungarian_algorithm(SWAP, energies), [1, 0])


def test_energy_aware_cost_scaling_matches_cpp_branches():
    overlap = np.ones((2, 2), dtype=complex)
    energies = np.diag([0.0, 2.0])
    np.testing.assert_allclose(make_cost_matrix(overlap, energies), np.ones((2, 2)))
    gaussian = make_cost_matrix(overlap, energies, alpha=0.5, scaling_function=1)
    symmetric = make_cost_matrix(overlap, energies, alpha=0.5, scaling_function=2)
    uphill = make_cost_matrix(overlap, energies, alpha=0.5, scaling_function=3)
    assert gaussian[0, 1] == pytest.approx(np.exp(-1.0))
    assert symmetric[0, 1] == pytest.approx(np.exp(-1.0))
    assert uphill[0, 1] == pytest.approx(1.0)
    assert uphill[1, 0] == pytest.approx(np.exp(-1.0))


def test_stochastic_complete_permutations_and_nonconvergence_behavior():
    np.testing.assert_array_equal(
        get_stochastic_reordering2(SWAP, np.random.default_rng(3)), [1, 0]
    )
    impossible = np.ones((2, 2), dtype=complex)
    identity = get_stochastic_reordering3(
        impossible, np.random.default_rng(1), convergence=0, max_number_of_attempts=0
    )
    np.testing.assert_array_equal(identity, [0, 1])
    with pytest.raises(RuntimeError):
        get_stochastic_reordering3(
            impossible, np.random.default_rng(1), convergence=1,
            max_number_of_attempts=0,
        )


def test_permutation_matrix_and_active_state_mapping_use_cpp_convention():
    matrix = permutation_matrix([2, 0, 1])
    np.testing.assert_allclose(matrix[[2, 0, 1], [0, 1, 2]], 1.0)
    np.testing.assert_array_equal(permute_states([[1, 0], [0, 1]], [0, 1]), [1, 1])


def test_force_cost_detects_predicted_gap_sign_change():
    forces = np.array([[0.0], [2.0]])
    previous_energies = np.diag([0.0, 0.1])
    cost = compute_force_cost_matrix(
        forces, forces, previous_energies, previous_energies,
        momentum=[0.0], inverse_mass=[1.0], dt=1.0,
    )
    assert cost[0, 1] == pytest.approx(1.0)
    assert cost[1, 0] == pytest.approx(1.0)


def test_projector_combines_reordering_and_phase_correction():
    overlap = np.array([[0.0, 1.0j], [-1.0, 0.0]])
    params = DynControlParams(
        state_tracking_algo=2,
        do_phase_correction=1,
        phase_correction_tol=1.0e-8,
    )
    projector = compute_projector(params, np.diag([0.0, 0.1]), overlap)
    corrected = overlap @ projector
    np.testing.assert_allclose(np.diag(corrected), [1.0, 1.0])
    np.testing.assert_allclose(projector.conj().T @ projector, np.eye(2))


def test_ld_and_svd_projection_updates_are_unitary():
    overlap = np.array([[0.9, 0.2], [-0.1, 0.8]], dtype=complex)
    for algorithm in (-1, 5, 6):
        projector = update_projection(
            {"state_tracking_algo": algorithm}, overlap,
            energies_current=np.diag([0.0, 0.1]),
        )
        np.testing.assert_allclose(projector.conj().T @ projector, np.eye(2), atol=1e-12)

