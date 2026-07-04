import numpy as np

from libra_py.dyn.models import GLVCModel

from _test_utils import assert_engine_builds_adiabatic


PARAMS = {
    "nstates": 2,
    "num_osc": 2,
    "Ham": [[0.0, 0.002], [0.002, 0.01]],
    "omega": [[0.01, 0.02], [0.015, 0.025]],
    "coupl": [[0.001, -0.002], [-0.001, 0.003]],
    "mass": [1.0, 2.0],
    "coupling_scaling": [1.0, -1.0],
}


def test_glvc_matches_reference_values():
    Q = np.asarray([[0.2], [-0.3]])
    result = GLVCModel(params=PARAMS).evaluate(Q)
    expected_00 = 0.5 * 1.0 * 0.01**2 * 0.2**2 + 0.001 * 0.2 + 0.5 * 2.0 * 0.02**2 * (-0.3) ** 2 + (-0.002) * (-0.3)
    expected_11 = 0.01 + 0.5 * 1.0 * 0.015**2 * 0.2**2 + (-0.001) * 0.2 * -1.0 + 0.5 * 2.0 * 0.025**2 * (-0.3) ** 2 + 0.003 * (-0.3) * -1.0
    assert np.allclose(result["H_dia"][0, 0, 0], expected_00)
    assert np.allclose(result["H_dia"][0, 1, 1], expected_11)
    assert np.allclose(result["H_dia"][0, 0, 1], 0.002)


def test_glvc_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(GLVCModel(params=PARAMS), np.asarray([[0.2, 0.3], [-0.3, 0.1]]))
