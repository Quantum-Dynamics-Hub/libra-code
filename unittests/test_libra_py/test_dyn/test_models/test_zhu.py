import numpy as np

from libra_py.dyn.models import ZhuDualLZSModel, ZhuDualRZDModel

from _test_utils import assert_engine_builds_adiabatic


def test_zhu_dual_lzs_matches_reference_values():
    x = np.array([-1.0, 2.0])
    result = ZhuDualLZSModel().evaluate(np.asarray([x]))
    assert np.allclose(result["H_dia"][:, 1, 1], 0.03 - 0.1 * np.exp(-0.28 * x * x))
    assert np.allclose(result["H_dia"][:, 0, 1], 0.01 * np.exp(-0.06 * x * x))


def test_zhu_variants_hamiltonian_engine_path():
    q = np.asarray([[-1.0, 2.0]])
    assert_engine_builds_adiabatic(ZhuDualLZSModel(), q)
    assert_engine_builds_adiabatic(ZhuDualRZDModel(), q)
