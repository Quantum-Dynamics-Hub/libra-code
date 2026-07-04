import numpy as np

from libra_py.dyn.models import FerrettiModel

from _test_utils import assert_engine_builds_adiabatic


def test_ferretti_matches_reference_values():
    Q = np.asarray([[0.5], [-0.2]])
    result = FerrettiModel().evaluate(Q)
    h12 = 0.005 * -0.2 * np.exp(-3.0 * (0.5 - 3.0) ** 2) * np.exp(-1.5 * (-0.2) ** 2)
    assert np.allclose(result["H_dia"][0, 0, 1], h12)


def test_ferretti_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(FerrettiModel(), np.asarray([[0.5, 1.0], [-0.2, 0.1]]))
