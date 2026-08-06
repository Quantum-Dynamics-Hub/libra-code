import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.hopping.hop_acceptance import (
    Boltz_cl_prob_up,
    Boltz_quant_prob,
    HO_prob,
    HO_prob_up,
    accept_hops,
    boltz_factor,
    can_rescale_along_vector,
    handle_hops_nuclear,
    rescale_along_vector,
    where_can_we_hop,
)


def test_rescaling_discriminant_detects_allowed_and_frustrated_hops():
    assert can_rescale_along_vector(1.0, 0.0, [0.0], [1.0], [1.0])
    assert not can_rescale_along_vector(0.0, 1.0, [1.0], [1.0], [1.0])


def test_rescale_along_vector_conserves_energy():
    momentum = np.array([0.0])
    rescale_along_vector(1.0, 0.0, momentum, [1.0], [1.0])
    np.testing.assert_allclose(0.5 * momentum[0] ** 2, 1.0)


def test_frustrated_rescaling_can_reverse_direction_component():
    momentum = np.array([1.0])
    rescale_along_vector(0.0, 1.0, momentum, [1.0], [1.0], do_reverse=True)
    np.testing.assert_allclose(momentum, [-1.0])


def test_quantum_and_classical_boltzmann_probabilities():
    probabilities = Boltz_quant_prob([0.0, 0.01], 300.0)
    np.testing.assert_allclose(probabilities.sum(), 1.0)
    assert probabilities[0] > probabilities[1] > 0.0
    assert 0.0 < Boltz_cl_prob_up(0.01, 300.0) < 1.0
    assert boltz_factor(0.0, 0.01, 300.0, 1) == 1.0
    np.testing.assert_allclose(
        boltz_factor(0.01, 0.0, 300.0, 3), probabilities[1]
    )


def test_harmonic_oscillator_probabilities_and_output_vector():
    output = []
    total, individual = HO_prob([0.01, 0.02], [0, 1], 300.0, output)
    np.testing.assert_allclose(total, np.prod(individual))
    np.testing.assert_allclose(output, individual)

    total_up, individual_up = HO_prob_up([0.01, 0.02], [4, 7], 300.0)
    np.testing.assert_allclose(total_up, np.prod(individual_up))


def test_accept_all_and_energy_conservation_options():
    energies = np.array([0.0, 1.0])
    proposed = [1, 0]
    initial = [0, 1]
    np.testing.assert_array_equal(
        accept_hops(DynControlParams(hop_acceptance_algo=0), proposed, initial, energies),
        proposed,
    )
    result = accept_hops(
        DynControlParams(hop_acceptance_algo=10),
        proposed,
        initial,
        energies,
        momenta=[[1.0], [0.0]],
        inverse_mass=[1.0],
    )
    np.testing.assert_array_equal(result, [0, 0])


def test_derivative_coupling_acceptance_uses_rescaling_feasibility():
    params = DynControlParams(hop_acceptance_algo=20)
    dc = np.array([[[[0.0, 1.0], [-1.0, 0.0]]]])
    rejected = accept_hops(
        params,
        [1],
        [0],
        [0.0, 1.0],
        momenta=[[1.0]],
        inverse_mass=[1.0],
        dc1_adi=dc,
    )
    accepted = accept_hops(
        params,
        [1],
        [0],
        [0.0, 1.0],
        momenta=[[2.0]],
        inverse_mass=[1.0],
        dc1_adi=dc,
    )
    np.testing.assert_array_equal(rejected, [0])
    np.testing.assert_array_equal(accepted, [1])


def test_where_can_we_hop_enumerates_energy_allowed_targets():
    possible = where_can_we_hop(
        0,
        DynControlParams(hop_acceptance_algo=10),
        [0],
        [0.0, 0.2, 1.0],
        momenta=[[1.0]],
        inverse_mass=[1.0],
    )
    assert possible == [1]


def test_uniform_nuclear_rescaling_and_tcnbra_update():
    momentum = np.array([[2.0]])
    handle_hops_nuclear(
        DynControlParams(momenta_rescaling_algo=100),
        momentum,
        [1.0],
        [1],
        [0],
        [0.0, 1.0],
    )
    np.testing.assert_allclose(momentum, [[np.sqrt(2.0)]])

    kinetic = np.array([2.0])
    handle_hops_nuclear(
        DynControlParams(momenta_rescaling_algo=40),
        momentum,
        [1.0],
        [1],
        [0],
        [0.0, 0.5],
        tcnbra_ekin=kinetic,
    )
    np.testing.assert_allclose(kinetic, [1.5])
