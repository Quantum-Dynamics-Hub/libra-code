import numpy as np

from libra_py.dyn.models import BeswickJortnerModel
from libra_py.units import Angst, ev2Ha

from _test_utils import assert_engine_builds_adiabatic


def test_beswick_jortner_matches_reference_formula():
    model = BeswickJortnerModel()
    q = np.asarray(
        [
            [1.20 * Angst, 1.35 * Angst],
            [2.30 * Angst, 2.55 * Angst],
        ]
    )

    result = model.evaluate(q)
    expected_h, expected_dh = _reference(q)

    assert result["H_dia"].shape == (2, 1, 1)
    assert result["dH_dia"].shape == (2, 2, 1, 1)
    assert np.allclose(result["H_dia"], expected_h)
    assert np.allclose(result["dH_dia"], expected_dh)
    assert np.allclose(result["S_dia"], np.ones((2, 1, 1)))


def test_beswick_jortner_works_through_hamiltonian_engine():
    q = np.asarray(
        [
            [1.20 * Angst, 1.35 * Angst, 1.50 * Angst],
            [2.30 * Angst, 2.55 * Angst, 2.80 * Angst],
        ]
    )
    assert_engine_builds_adiabatic(BeswickJortnerModel(), q)


def _reference(q):
    K = 74.4434 * ev2Ha / (Angst * Angst)
    r0 = 1.2327 * Angst
    T0 = 4.36994 * ev2Ha
    A = 200000.0 * ev2Ha
    a = 6.68 / Angst
    m_c = 12.0
    m_n = 14.0
    fraction = m_c / (m_c + m_n)
    r, R = q[0], q[1]
    exp_term = np.exp(a * (fraction * r - R))
    repulsive = A * exp_term
    value = T0 + 0.5 * K * (r - r0) ** 2 + repulsive
    d_dr = K * (r - r0) + a * fraction * repulsive
    d_dR = -a * repulsive
    h = np.zeros((q.shape[1], 1, 1), dtype=complex)
    dh = np.zeros((q.shape[1], 2, 1, 1), dtype=complex)
    h[:, 0, 0] = value
    dh[:, 0, 0, 0] = d_dr
    dh[:, 1, 0, 0] = d_dR
    return h, dh
