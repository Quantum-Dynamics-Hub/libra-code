import numpy as np

from libra_py.dyn.core.manager import TrajectoryManager
from libra_py.dyn.core.storage import TensorStorage


def test_manager_create_spawn_copy_and_lifecycle():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=5,
    )
    manager = TrajectoryManager(storage)

    traj = manager.create_trajectory(metadata={"label": "demo"})
    parent = manager.create_tbf(traj, idx=0, metadata={"role": "parent"})
    parent.view.q[:] = [0.0, 0.5]
    parent.view.ampl_adi[:] = [1.0, 0.0, 0.0]

    child = manager.spawn_tbf(
        traj,
        parent=parent,
        spawn_time=2.5,
        metadata={"role": "spawned"},
    )

    assert child.parent_id == parent.id
    assert child.spawn_time == 2.5
    assert np.allclose(child.view.q, parent.view.q)
    assert np.allclose(child.view.ampl_adi, parent.view.ampl_adi)
    assert manager.active_tbf_ids(traj) == [0, 1]
    assert [tbf.id for tbf in manager.active_tbfs(traj)] == [0, 1]

    manager.deactivate_tbf(child)
    assert manager.active_tbf_ids(traj) == [0]
    assert not bool(storage.alive[traj.id, child.tbf_id])

    manager.activate_tbf(child.id)
    assert manager.active_tbf_ids(traj) == [0, 1]

    removed = manager.remove_tbf(child.id)
    assert removed.id == child.id
    assert child.id not in manager.tbfs
    assert traj.tbf_ids == [0]


def test_manager_trajectory_lifecycle_and_compaction():
    storage = TensorStorage(
        backend=np,
        ntraj=3,
        ndof=2,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=5,
    )
    manager = TrajectoryManager(storage)

    traj0 = manager.create_trajectory()
    traj1 = manager.create_trajectory()
    tbf0 = manager.create_tbf(traj0, idx=0)
    tbf1 = manager.create_tbf(traj1, idx=0)
    spawned = manager.spawn_tbf(traj0, parent=tbf0)

    manager.deactivate_trajectory(traj1, deactivate_tbfs=True)
    assert not traj1.active
    assert manager.active_tbf_ids(traj1) == []

    manager.activate_trajectory(traj1, activate_tbfs=True)
    assert traj1.active
    assert manager.active_tbf_ids(traj1) == [0]

    manager.remove_tbf(spawned)
    slot_map = manager.compact_storage()

    assert slot_map == {0: 0}
    assert storage.ntbf_capacity == 1
    assert tbf0.tbf_id == 0
    assert tbf1.tbf_id == 0

    summary = manager.summary()
    assert summary["n_trajectories"] == 2
    assert summary["n_tbfs"] == 2
