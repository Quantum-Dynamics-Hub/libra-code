"""Optional Gaussian fields follow the existing TBF storage lifecycle."""
import numpy as np
import pytest
from libra_py.dyn.backends import NumpyBackend, TorchBackend
from libra_py.dyn.backends.arrays import to_numpy
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.manager import TrajectoryManager


@pytest.mark.parametrize('backend', ['numpy', 'torch'])
def test_gaussian_allocation_spawn_and_compact(backend):
    if backend == 'torch':
        pytest.importorskip('torch')
    storage = TensorStorage(NumpyBackend() if backend == 'numpy' else TorchBackend(),
                            ntraj=1, ndof=2, nstates=2, ntbf_initial=1, ntbf_capacity=1)
    assert storage.gaussian_alpha is None
    manager = TrajectoryManager(storage)
    traj = manager.create_trajectory()
    parent = manager.create_tbf(traj, 0)
    with pytest.raises(AttributeError, match='not been allocated'):
        _ = parent.view.gaussian_alpha
    storage.allocate_gaussian_vars()
    parent.view.gaussian_alpha[:] = storage.backend.asarray([2., 3.])
    parent.view.gaussian_phase[:] = storage.backend.asarray([.4, .6])
    storage.allocate_gaussian_vars()  # allocation is idempotent
    child = manager.spawn_tbf(traj, parent=parent)
    assert storage.ntbf_capacity == 2
    np.testing.assert_array_equal(to_numpy(child.view.gaussian_alpha), [2., 3.])
    np.testing.assert_array_equal(to_numpy(child.view.gaussian_phase), [.4, .6])
    manager.deactivate_tbf(parent)
    manager.compact_storage()
    assert child.tbf_id == 0
    np.testing.assert_array_equal(to_numpy(child.view.gaussian_alpha), [2., 3.])
