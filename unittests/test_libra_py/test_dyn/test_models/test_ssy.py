import numpy as np

from libra_py.dyn.models import SSYModel

from _test_utils import assert_engine_builds_adiabatic


def test_ssy_matches_reference_values():
    Q = np.asarray([[0.5], [-0.2]])
    result = SSYModel().evaluate(Q)
    z = 0.015 * np.exp(-0.06 * (0.25 * (0.5 - 0.2) ** 2 + 0.75 * (0.5 + 0.2) ** 2))
    assert np.allclose(result["H_dia"][0, 0, 1], z)


def test_ssy_hamiltonian_engine_path():
    assert_engine_builds_adiabatic(SSYModel(), np.asarray([[0.5, 0.0], [-0.2, 0.1]]))
