import numpy as np

from libra_py.dyn.models import LibraModel1

from _test_utils import assert_engine_builds_adiabatic


def test_libra_model1_matches_reference_values():
    x = np.array([-1.0, 0.5, 2.0])
    params = {"x0": 1.5, "k": 0.02, "D": -0.01, "V": 0.003}
    result = LibraModel1(params=params).evaluate(np.asarray([x]))
    expected = np.zeros((3, 2, 2), dtype=complex)
    expected[:, 0, 0] = params["k"] * x * x
    expected[:, 1, 1] = params["k"] * (x - params["x0"]) ** 2 + params["D"]
    expected[:, 0, 1] = params["V"]
    expected[:, 1, 0] = params["V"]

    assert np.allclose(result["H_dia"], expected)
    assert np.allclose(result["dH_dia"][:, 0, 0, 0], 2.0 * params["k"] * x)
    assert np.allclose(result["dH_dia"][:, 0, 1, 1], 2.0 * params["k"] * (x - params["x0"]))


def test_libra_model1_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(LibraModel1(), np.asarray([[-1.0, 1.0]]))
