import numpy as np

from libra_py.dyn.models import Holstein2Model, Holstein3Model, Holstein4Model, Holstein5Model

from _test_utils import assert_engine_builds_adiabatic


BASE = {"E_n": [0.0, 0.2, 0.3], "x_n": [0.0, 1.0, 2.0], "k_n": [0.1, 0.2, 0.3]}


def test_holstein2_matches_reference_values():
    x = np.array([-0.5, 1.5])
    result = Holstein2Model(params={**BASE, "V": 0.01}).evaluate(np.asarray([x]))
    assert np.allclose(result["H_dia"][:, 0, 0], 0.5 * BASE["k_n"][0] * x * x)
    assert np.allclose(result["H_dia"][:, 0, 1], 0.01)


def test_holstein_variants_hamiltonian_engine_path():
    q = np.asarray([[-0.5, 1.5]])
    assert_engine_builds_adiabatic(Holstein2Model(params={**BASE, "V": 0.01}), q)
    assert_engine_builds_adiabatic(Holstein3Model(params={**BASE, "V_n": [0.01, 0.002]}), q)
    assert_engine_builds_adiabatic(Holstein4Model(params={**BASE, "V": [[0.0, 0.01, 0.002], [0.01, 0.0, 0.01], [0.002, 0.01, 0.0]]}), q)
    assert_engine_builds_adiabatic(Holstein5Model(params={**BASE, "V": [[0.0, 0.01, 0.002], [0.01, 0.0, 0.01], [0.002, 0.01, 0.0]], "alpha": [[0.0, 0.2, 0.1], [0.2, 0.0, 0.2], [0.1, 0.2, 0.0]], "x_nm": [[0.0, 0.5, 1.0], [0.5, 0.0, 1.5], [1.0, 1.5, 0.0]]}), q)
