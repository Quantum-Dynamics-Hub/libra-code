import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.hopping.hop_proposal_fssh3 import (
    adjust_signs,
    find_best_matrix,
    hopping_probabilities_fssh3,
)
from libra_py.dyn.hopping.hop_proposal import (
    TSH_METHODS,
    hop,
    hop_proposal_probabilities,
    hopping_probabilities_fssh,
    hopping_probabilities_fssh2,
    hopping_probabilities_gfsh_orig,
    hopping_probabilities_lz,
    hopping_probabilities_mash,
    hopping_probabilities_mssh,
    hopping_probabilities_zn,
    propose_hops,
)


def test_fssh_matches_cpp_reference_case():
    params = DynControlParams(dt=41.0, Temperature=300.0, use_boltz_factor=1)
    hvib = np.array([[0.0, -0.1j], [0.1j, 0.01]], dtype=complex)
    density = np.array([[0.5, 1.0], [0.0, 0.5]], dtype=complex)

    result = hopping_probabilities_fssh(params, density, hvib, 0)

    np.testing.assert_allclose(result, [0.9997799598333882, 0.00022004016661186391])


def test_fssh_amplitudes_return_transition_matrix():
    params = DynControlParams(dt=0.1)
    hvib = np.array([[0.0, -0.1j], [0.1j, 0.0]], dtype=complex)
    result = hopping_probabilities_fssh(params, [2**-0.5, 2**-0.5], hvib)

    assert result.shape == (2, 2)
    np.testing.assert_allclose(result.sum(axis=1), 1.0)


def test_population_difference_methods():
    old = np.diag([0.8, 0.2]).astype(complex)
    new = np.diag([0.6, 0.4]).astype(complex)

    np.testing.assert_allclose(
        hopping_probabilities_gfsh_orig({}, new, old, 0), [2.0 / 3.0, 1.0 / 3.0]
    )
    np.testing.assert_allclose(
        hopping_probabilities_fssh2({"fssh2_revision": 0}, new, old, 0),
        [0.75, 0.25],
    )


def test_mssh_boltzmann_scales_uphill_hop():
    params = DynControlParams(Temperature=300.0, use_boltz_factor=1)
    density = np.diag([0.5, 0.5]).astype(complex)
    hvib = np.diag([0.0, 0.01]).astype(complex)

    result = hopping_probabilities_mssh(params, density, hvib, 0)

    assert 0.0 < result[1] < 0.5
    np.testing.assert_allclose(result.sum(), 1.0)


def test_mash_selects_largest_population():
    np.testing.assert_array_equal(
        hopping_probabilities_mash({}, np.diag([0.2, 0.7, 0.1])), [0.0, 1.0, 0.0]
    )


def test_hop_and_batch_proposal_are_deterministic_with_seed():
    assert hop(0, [0.2, 0.8], 0.1) == 0
    assert hop(0, [0.2, 0.8], 0.5) == 1
    result = propose_hops([[0.2, 0.8], [0.9, 0.1]], [0, 1], np.random.default_rng(2))
    np.testing.assert_array_equal(result, [1, 0])


def test_dispatcher_uses_control_parameter_method():
    params = DynControlParams(tsh_method=6)
    density = np.diag([0.3, 0.7]).astype(complex)
    result = hop_proposal_probabilities(params, density, np.eye(2), 0)

    np.testing.assert_array_equal(result, [0.0, 1.0])


def test_lz_diabatic_probability_at_gap_crossing():
    current = {
        "ham_dia": np.array([[0.1, 0.01], [0.01, -0.1]]),
        "d1ham_dia": np.array([[[1.0, 0.0], [0.0, -1.0]]]),
    }
    previous = {"ham_dia": np.diag([-0.1, 0.1])}

    result = hopping_probabilities_lz(current, previous, 0, 0, [1.0], [1.0])

    expected_hop = np.exp(-2.0 * np.pi * 0.01**2 / 2.0)
    np.testing.assert_allclose(result, [1.0 - expected_hop, expected_hop])


def test_lz_adiabatic_probability_and_nac_crossing_detection():
    current = {
        "ham_dia": np.diag([0.1, -0.1]),
        "ham_adi": np.diag([0.0, 0.02]),
        "nac_adi": np.array([[0.0, 0.1], [-0.1, 0.0]]),
    }
    previous = {
        "ham_dia": np.diag([-0.1, 0.1]),
        "nac_adi": np.array([[0.0, -0.1], [0.1, 0.0]]),
    }

    result = hopping_probabilities_lz(current, previous, 0, 2, [1.0], [1.0])

    expected_hop = np.exp(-0.25 * np.pi * 0.02 / 0.1)
    np.testing.assert_allclose(result, [1.0 - expected_hop, expected_hop])


def test_lz_requires_a_crossing():
    current = {
        "ham_dia": np.array([[0.2, 0.01], [0.01, 0.0]]),
        "d1ham_dia": np.array([[[1.0, 0.0], [0.0, -1.0]]]),
    }
    previous = {"ham_dia": np.diag([0.1, 0.0])}

    np.testing.assert_array_equal(
        hopping_probabilities_lz(current, previous, 0, 0, [1.0], [1.0]),
        [1.0, 0.0],
    )


def test_zn_probability_matches_multidimensional_formula():
    current = {
        "ham_dia": np.array([[0.1, 0.1], [0.1, -0.1]]),
        "forces_adi": np.array([[1.0, -1.0]]),
    }
    previous = {"ham_dia": np.diag([-0.1, 0.1])}

    result = hopping_probabilities_zn(current, previous, 0, 0, [0.0], [1.0])

    a2 = 0.0625 * 1.0 * 2.0 / 0.1**3
    b2 = 0.5 * 2.0 / (1.0 * 0.1)
    expected_hop = np.exp(
        -(0.25 * np.pi / np.sqrt(a2))
        * np.sqrt(2.0 / (b2 + np.sqrt(b2**2 - 1.0)))
    )
    np.testing.assert_allclose(result, [1.0 - expected_hop, expected_hop])


def test_dispatcher_calls_lz_with_hamiltonian_history():
    params = DynControlParams(tsh_method=3, rep_lz=0)
    current = {
        "ham_dia": np.array([[0.1, 0.01], [0.01, -0.1]]),
        "d1ham_dia": np.array([[[1.0, 0.0], [0.0, -1.0]]]),
    }
    previous = {"ham_dia": np.diag([-0.1, 0.1])}

    direct = hopping_probabilities_lz(current, previous, 0, 0, [1.0], [1.0])
    dispatched = hop_proposal_probabilities(
        params,
        None,
        None,
        0,
        ham=current,
        ham_prev=previous,
        momentum=[1.0],
        inverse_mass=[1.0],
    )

    np.testing.assert_allclose(dispatched, direct)


def test_fssh3_uses_observed_active_population_loss():
    params = DynControlParams(
        dt=1.0,
        fssh3_approach_option=0,
        fssh3_dt=0.01,
        fssh3_max_steps=2000,
        fssh3_err_tol=1.0e-12,
    )
    old = np.diag([0.8, 0.2]).astype(complex)
    new = np.diag([0.6, 0.4]).astype(complex)
    errors = [0.0] * 5

    result = hopping_probabilities_fssh3(params, new, old, 0, errors)

    np.testing.assert_allclose(result, [0.75, 0.25])
    assert errors[0] == errors[1]
    assert errors[0] < 1.0e-12
    np.testing.assert_array_equal(errors[2:], [0.0, 0.0, 0.0])


def test_fssh3_zero_outflux_stays_on_active_state():
    density = np.diag([0.7, 0.3]).astype(complex)
    result = hopping_probabilities_fssh3({}, density, density, 0)
    np.testing.assert_array_equal(result, [1.0, 0.0])


def test_fssh3_sign_helpers_return_best_candidate():
    matrix = np.array([[1.0, 2.0], [-3.0, 4.0]])
    adjusted = adjust_signs(matrix)
    np.testing.assert_array_equal(adjusted, [[1.0, 2.0], [-3.0, 4.0]])

    old = np.array([1.0, 0.0])
    new = np.array([-1.0, 3.0])
    best = find_best_matrix(new, old, matrix)
    np.testing.assert_allclose(best @ old, new)


def test_dispatcher_calls_fssh3():
    params = DynControlParams(
        tsh_method=8,
        dt=1.0,
        fssh3_dt=0.01,
        fssh3_max_steps=2000,
        fssh3_err_tol=1.0e-12,
    )
    old = np.diag([0.8, 0.2]).astype(complex)
    new = np.diag([0.6, 0.4]).astype(complex)

    result = hop_proposal_probabilities(params, new, np.eye(2), 0, old)

    np.testing.assert_allclose(result, [0.75, 0.25])


def test_cpp_method_numbering_is_complete():
    assert tuple(TSH_METHODS) == (-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)


def test_no_hop_option_returns_active_state():
    params = DynControlParams(tsh_method=-1)
    result = hop_proposal_probabilities(params, np.eye(2), np.eye(2), 1)
    np.testing.assert_array_equal(result, [0.0, 1.0])


def test_dish_option_is_an_explicit_placeholder():
    params = DynControlParams(tsh_method=5)
    with np.testing.assert_raises_regex(NotImplementedError, "event scheduler"):
        hop_proposal_probabilities(params, np.eye(2), np.eye(2), 0)


def test_dispatcher_batches_lz_probabilities():
    params = DynControlParams(tsh_method=3, rep_lz=0)
    current = {
        "ham_dia": np.array([[0.1, 0.01], [0.01, -0.1]]),
        "d1ham_dia": np.array([[[1.0, 0.0], [0.0, -1.0]]]),
    }
    previous = {"ham_dia": np.diag([-0.1, 0.1])}

    result = hop_proposal_probabilities(
        params,
        None,
        None,
        [0, 0],
        ham=[current, current],
        ham_prev=[previous, previous],
        momentum=np.array([[1.0, 2.0]]),
        inverse_mass=[1.0],
    )

    assert result.shape == (2, 2)
    assert result[1, 1] > result[0, 1]
