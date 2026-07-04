import numpy as np

from libra_py.dyn.models import MartensModel1, MartensModel2

from _test_utils import assert_engine_builds_adiabatic


def test_martens_model1_matches_reference_value():
    Q = np.asarray([[0.2], [-0.4]])
    result = MartensModel1().evaluate(Q)
    expected = 0.00625 / np.cosh(0.4) ** 2 + 0.5 * 0.0106 * (-0.4) ** 2
    assert np.allclose(result["H_dia"][0, 0, 0], expected)


def test_martens_model2_matches_reference_value():
    Q = np.asarray([[0.2], [-0.4]])
    result = MartensModel2().evaluate(Q)
    shift = -0.4 + 0.4 * (0.2 * 0.2 - 1.0)
    expected = 0.00625 / np.cosh(0.4) ** 2 + 0.5 * 0.0106 * shift * shift
    assert np.allclose(result["H_dia"][0, 0, 0], expected)


def test_martens_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(MartensModel2(), np.asarray([[0.2, 0.3], [-0.4, 0.1]]))
