import numpy as np

from libra_py.dyn.models import LVCModel

from _test_utils import assert_engine_builds_adiabatic


PARAMS = {"Delta1": 0.0, "Delta2": 0.1, "omega": [0.01, 0.02], "d1": [0.0, 0.0], "d2": [0.01, -0.01], "coup": [0.001, 0.002], "mass": [1.0, 4.0]}


def test_lvc_matches_reference_values():
    Q = np.asarray([[0.2], [-0.4]])
    result = LVCModel(params=PARAMS).evaluate(Q)
    expected_coupling = np.sqrt(1.0) * 0.001 * 0.2 + np.sqrt(4.0) * 0.002 * -0.4
    assert np.allclose(result["H_dia"][0, 0, 1], expected_coupling)
    assert np.allclose(result["dH_dia"][0, 1, 1, 1], 4.0 * 0.02 * 0.02 * -0.4 + np.sqrt(4.0) * -0.01)


def test_lvc_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(LVCModel(params=PARAMS), np.asarray([[0.2, 0.3], [-0.4, 0.1]]))
