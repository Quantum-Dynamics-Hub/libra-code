import numpy as np

from libra_py.dyn.decoherence import (
    coherence_intervals, collapse, decoherence_event, edc_rates,
    gu_franco, instantaneous_decoherence, project_out, schwartz_2, sdm,
)


def test_edc_rates_two_state_exact_value():
    h = np.diag([0.0, 0.2]).astype(complex)
    rates = edc_rates(h, 0.5, C_param=1.0, eps_param=0.1)
    assert np.allclose(rates, [[0.0, 1.0 / 6.0], [1.0 / 6.0, 0.0]])


def test_coherence_intervals_use_other_state_populations():
    c = np.sqrt([0.75, 0.25]).astype(complex)
    rates = np.array([[0.0, 2.0], [4.0, 0.0]])
    assert np.allclose(coherence_intervals(c, rates), [2.0, 1.0 / 3.0])


def test_sdm_decays_inactive_and_preserves_norm_and_phase():
    c = np.array([[np.sqrt(0.6), 1j * np.sqrt(0.4)]])
    rates = np.array([[0.0, 1.0], [1.0, 0.0]])
    result = sdm(c, 0.5, [0], rates)
    assert np.isclose(abs(result[0, 1]), np.sqrt(0.4) * np.exp(-0.5))
    assert np.isclose(np.vdot(result[0], result[0]).real, 1.0)
    assert np.isclose(result[0, 1].real, 0.0)


def test_projection_collapse_and_id_failed_hop():
    c = np.array([[np.sqrt(0.8), 1j*np.sqrt(0.2)]], complex)
    instantaneous_decoherence(c, [0], [1], [0], variant=3)
    assert np.allclose(c, [[1.0, 0.0]])
    c = np.array([[1j/np.sqrt(2), 1/np.sqrt(2)]], complex)
    collapse(c, 0, 0, collapse_option=0)
    assert np.allclose(c, [[1j, 0]])
    project_out(c, 0, 0)
    assert np.allclose(c, 0.0)


def test_schwartz_2_pair_rate_and_gu_franco_symmetry():
    forces = np.array([[[1.0, 0.0], [-1.0, 0.0], [0.0, 2.0]]])
    rates = schwartz_2(forces, [1.0, 4.0])
    assert np.isclose(rates[0, 0, 1], 1.0)
    assert np.allclose(rates, rates.swapaxes(1, 2))
    gf = gu_franco([[1/np.sqrt(2), 1j/np.sqrt(2)]], 0.1, 300.0)
    assert np.isclose(gf[0, 0, 1], gf[0, 1, 0])


def test_decoherence_event_selects_only_expired_state():
    events = decoherence_event([[0.2, 2.0, 0.1]], [[1.0, 1.0, 1.0]], option=0,
                               rng=np.random.default_rng(2))
    assert events.tolist() == [1]
