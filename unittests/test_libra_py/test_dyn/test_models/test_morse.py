import numpy as np

from libra_py.dyn.models import MorseModel

from _test_utils import assert_engine_builds_adiabatic


PARAMS = {"E": [0.0, 0.1], "D": [0.2, 0.3], "alpha": [0.4, 0.5], "x_n": [0.0, 1.0], "V": [[0.0, 0.01], [0.01, 0.0]], "beta": [[0.0, 0.2], [0.2, 0.0]], "x_nm": [[0.0, 0.25], [0.25, 0.0]]}


def test_morse_matches_reference_values():
    x = np.array([-0.5, 1.5])
    result = MorseModel(params=PARAMS).evaluate(np.asarray([x]))
    bo = np.exp(-0.4 * x)
    assert np.allclose(result["H_dia"][:, 0, 0], 0.2 * (1.0 - bo) ** 2)
    assert np.allclose(result["H_dia"][:, 0, 1], 0.01 * np.exp(-0.2 * (x - 0.25) ** 2))


def test_morse_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(MorseModel(params=PARAMS), np.asarray([[-0.5, 1.5]]))
