import numpy as np
import pytest

from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.views import TensorView


def test_tensor_view_writable_core_slices():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=3,
        nstates=2,
        ntbf_initial=1,
    )
    view = TensorView(storage, traj_id=1, tbf_id=0)

    view.q[:] = [0.1, 0.2, 0.3]
    view.ampl_adi[:] = [1.0 + 0.0j, 0.0 + 0.0j]
    view.ham_adi[:] = np.diag([0.0, 0.25])
    view.set("act_states", 1)

    assert np.allclose(storage.q[1, 0], [0.1, 0.2, 0.3])
    assert storage.ampl_adi[1, 0, 0] == 1.0
    assert storage.ham_adi[1, 0, 1, 1] == 0.25
    assert storage.act_states[1, 0] == 1


def test_tensor_view_grouped_helpers_and_optional_fields():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=3,
        nstates=2,
        ntbf_initial=1,
    )
    view = TensorView(storage, traj_id=1, tbf_id=0)

    assert view.maybe("dR") is None
    with pytest.raises(AttributeError):
        _ = view.dR

    nuclear = view.nuclear()
    electronic = view.electronic()
    hamiltonian = view.hamiltonian()

    nuclear["q"][:] = [1.0, 2.0, 3.0]
    electronic["q_mm"][:] = [0.5, 0.6]
    hamiltonian["cum_phase_corr"][:] = np.eye(2)

    assert np.allclose(view.q, [1.0, 2.0, 3.0])
    assert np.allclose(view.q_mm, [0.5, 0.6])
    assert np.allclose(view.cum_phase_corr, np.eye(2))

    storage.allocate_afssh()
    view.dR[0, 0, 1] = 0.123
    assert storage.dR[1, 0, 0, 0, 1] == 0.123
    assert "dR" in view.optional()
