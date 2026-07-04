import numpy as np

from libra_py.dyn.models import FaistLevineModel, faist_levine_lii_params, faist_levine_nai_params

from _test_utils import assert_engine_builds_adiabatic


def test_faist_levine_matches_reference_coupling():
    r = np.asarray([5.0])
    params = faist_levine_nai_params()
    result = FaistLevineModel(params=params).evaluate(np.asarray([r]))
    h01 = params["A"] * np.exp(-r / params["rho"])
    assert np.allclose(result["H_dia"][:, 0, 1], h01)
    assert np.allclose(result["dH_dia"][:, 0, 0, 1], -h01 / params["rho"])


def test_faist_levine_variants_hamiltonian_engine_path():
    q = np.asarray([[5.0, 6.0]])
    assert_engine_builds_adiabatic(FaistLevineModel(params=faist_levine_nai_params()), q)
    assert_engine_builds_adiabatic(FaistLevineModel(params=faist_levine_lii_params()), q)
