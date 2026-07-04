import numpy as np

from libra_py.dyn.models import HenonHeilesModel

from _test_utils import assert_engine_builds_adiabatic


def test_henon_heiles_matches_reference_value():
    Q = np.asarray([[0.2], [-0.4]])
    result = HenonHeilesModel(params={"lam": 0.15}).evaluate(Q)
    x, y, lam = 0.2, -0.4, 0.15
    r2 = x * x + y * y
    expected = 0.5 * r2 + lam * (x * y * y - x**3 / 3.0) + lam * lam * r2 * r2 / 16.0
    assert np.allclose(result["H_dia"][0, 0, 0], expected)


def test_henon_heiles_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(HenonHeilesModel(), np.asarray([[0.2, 0.3], [-0.4, 0.1]]))
