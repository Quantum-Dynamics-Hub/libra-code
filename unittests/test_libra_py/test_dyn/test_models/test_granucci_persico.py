import numpy as np

from libra_py.dyn.models import GranucciPersicoModel1, GranucciPersicoModel2

from _test_utils import assert_engine_builds_adiabatic


def test_granucci_persico_model1_matches_reference_formula():
    x = np.asarray([-2.0, 0.0, 5.0])
    result = GranucciPersicoModel1().evaluate(np.asarray([x]))
    expected_h, expected_dh = _reference_model1(x)

    assert np.allclose(result["H_dia"], expected_h)
    assert np.allclose(result["dH_dia"], expected_dh)
    assert np.allclose(result["DC1_dia"], np.zeros_like(expected_dh))


def test_granucci_persico_model2_matches_reference_formula():
    q = np.asarray(
        [
            [2.0, 3.5, 5.0],
            [-0.5, 0.0, 0.75],
        ]
    )
    result = GranucciPersicoModel2().evaluate(q)
    expected_h, expected_dh = _reference_model2(q)

    assert result["H_dia"].shape == (3, 2, 2)
    assert result["dH_dia"].shape == (3, 2, 2, 2)
    assert np.allclose(result["H_dia"], expected_h)
    assert np.allclose(result["dH_dia"], expected_dh)


def test_granucci_persico_models_work_through_hamiltonian_engine():
    assert_engine_builds_adiabatic(GranucciPersicoModel1(), np.asarray([[-2.0, 0.0, 2.0]]))
    assert_engine_builds_adiabatic(
        GranucciPersicoModel2(),
        np.asarray([[2.0, 3.5, 5.0], [-0.5, 0.0, 0.75]]),
    )


def _reference_model1(x):
    a1, a2 = 0.006, 0.5
    alp1, alp2 = 0.2, 0.5
    dE, b, beta, gamma, x_c = 0.03, 0.013, 0.2, 0.0002, 5.0
    e1 = a1 * np.exp(-alp1 * x)
    e2 = a2 * np.exp(-alp2 * x)
    e3 = b * np.exp(-beta * (x - x_c) ** 2)
    h00 = e1 + dE
    h11 = e2
    h01 = e3 + gamma * np.sin(x) ** 2
    dh00 = -alp1 * e1
    dh11 = -alp2 * e2
    dh01 = -2.0 * beta * (x - x_c) * e3 + 2.0 * gamma * np.sin(x) * np.cos(x)
    return _two_state(x, h00, h11, h01, dh00, dh11, dh01)


def _reference_model2(q):
    x, y = q[0], q[1]
    D1, D2 = 0.015, 0.11
    delta1, delta2 = 0.05, 0.11
    alp1, alp2 = 1.0, 0.674
    beta1, beta2 = 0.5, 1.5
    gamma, K = 2.2295e-2, 0.09
    x1, x2, x3 = 3.9, 3.0, 5.0
    e1 = np.exp(-alp1 * (x - x1))
    e2 = np.exp(-alp2 * (x - x2))
    e3 = np.exp(-beta1 * (x - x3) ** 2 - beta2 * y * y)
    h00 = D1 * (e1 * e1 - 2.0 * e1) + delta1 + 0.5 * K * y * y
    h11 = D2 * (e2 * e2 - 2.0 * e2) + delta2 + 0.5 * K * y * y
    h01 = gamma * y * e3
    H = np.zeros((q.shape[1], 2, 2), dtype=complex)
    dH = np.zeros((q.shape[1], 2, 2, 2), dtype=complex)
    H[:, 0, 0] = h00
    H[:, 1, 1] = h11
    H[:, 0, 1] = h01
    H[:, 1, 0] = h01

    de1_dx = -alp1 * e1
    de2_dx = -alp2 * e2
    dH[:, 0, 0, 0] = 2.0 * D1 * (e1 - 1.0) * de1_dx
    dH[:, 0, 1, 1] = 2.0 * D2 * (e2 - 1.0) * de2_dx
    dH[:, 0, 0, 1] = -2.0 * beta1 * (x - x3) * h01
    dH[:, 0, 1, 0] = dH[:, 0, 0, 1]
    dH[:, 1, 0, 0] = K * y
    dH[:, 1, 1, 1] = K * y
    dH[:, 1, 0, 1] = gamma * e3 - 2.0 * beta2 * y * h01
    dH[:, 1, 1, 0] = dH[:, 1, 0, 1]
    return H, dH


def _two_state(x, h00, h11, h01, dh00, dh11, dh01):
    H = np.zeros((len(x), 2, 2), dtype=complex)
    dH = np.zeros((len(x), 1, 2, 2), dtype=complex)
    H[:, 0, 0] = h00
    H[:, 1, 1] = h11
    H[:, 0, 1] = h01
    H[:, 1, 0] = h01
    dH[:, 0, 0, 0] = dh00
    dH[:, 0, 1, 1] = dh11
    dH[:, 0, 0, 1] = dh01
    dH[:, 0, 1, 0] = dh01
    return H, dH
