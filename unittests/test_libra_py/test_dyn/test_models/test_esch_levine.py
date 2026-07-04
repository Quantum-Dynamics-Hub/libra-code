import numpy as np

from libra_py.dyn.models import EschLevineJCP2020Model, EschLevineLinearModel

from _test_utils import assert_engine_builds_adiabatic


def test_esch_levine_jcp2020_matches_reference_values():
    x = np.array([-0.5, 1.5])
    model = EschLevineJCP2020Model(params={"nstates": 3, "w0": 0.2, "w1": 0.05, "delta": 0.01, "V": 0.004, "eps": 0.02, "i_crit": 2})
    result = model.evaluate(np.asarray([x]))
    assert np.allclose(result["H_dia"][:, 0, 0], -0.2 * x)
    assert np.allclose(result["H_dia"][:, 2, 2], 0.05 * x - 0.02 - 0.02)
    assert np.allclose(result["H_dia"][:, 0, 2], 0.004)


def test_esch_levine_variants_hamiltonian_engine_path():
    q = np.asarray([[-0.5, 1.5]])
    assert_engine_builds_adiabatic(EschLevineJCP2020Model(params={"nstates": 3}), q)
    assert_engine_builds_adiabatic(EschLevineLinearModel(params={"nstates": 2, "V": [[0.0, 0.005], [0.005, 0.0]], "w": [[-0.1, 0.0], [0.0, 0.1]]}), q)
