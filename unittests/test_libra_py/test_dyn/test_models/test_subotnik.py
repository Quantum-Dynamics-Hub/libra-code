import numpy as np

from libra_py.dyn.models import SubotnikDoubleArchModel, SubotnikDumbbellModel

from _test_utils import assert_engine_builds_adiabatic


def test_subotnik_dumbbell_matches_reference_diagonals():
    x = np.array([-1.0, 2.0])
    result = SubotnikDumbbellModel(params={"Z": 1.0}).evaluate(np.asarray([x]))
    assert np.allclose(result["H_dia"][:, 0, 0], 0.0006)
    assert np.allclose(result["H_dia"][:, 1, 1], -0.0006)


def test_subotnik_variants_hamiltonian_engine_path():
    q = np.asarray([[-1.0, 2.0]])
    assert_engine_builds_adiabatic(SubotnikDumbbellModel(params={"Z": 1.0}), q)
    assert_engine_builds_adiabatic(SubotnikDoubleArchModel(params={"Z": 1.0}), q)
