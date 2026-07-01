import numpy as np

from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.tbf import TrajectoryBasisFunction
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.core.views import TensorView


def make_tbf(storage, object_id, traj_id, slot):
    return TrajectoryBasisFunction(
        id=object_id,
        view=TensorView(storage, traj_id=traj_id, tbf_id=slot),
        trajectory_id=traj_id,
    )


def test_trajectory_tracks_storage_slots_and_object_ids():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=4,
    )
    traj = Trajectory(id=1)

    tbf0 = make_tbf(storage, object_id=10, traj_id=1, slot=0)
    traj.add_tbf(tbf0)

    slot1 = storage.spawn(traj.id)
    tbf1 = make_tbf(storage, object_id=11, traj_id=1, slot=slot1)
    traj.add_tbf(tbf1)

    assert traj.tbf_ids == [0, 1]
    assert traj.tbf_object_ids == [10, 11]
    assert traj.tbf_object_by_slot == {0: 10, 1: 11}
    assert list(traj) == [0, 1]
    assert len(traj) == 2
    assert traj.active_tbf_ids(storage) == [0, 1]

    tbf1.deactivate()
    assert traj.active_tbf_ids(storage) == [0]

    traj.remove_tbf(slot1)
    traj.set_norm(0.998)
    traj.deactivate()

    assert traj.tbf_ids == [0]
    assert traj.tbf_object_ids == [10]
    assert not traj.active
    assert traj.norm == 0.998
    assert traj.to_dict()["metadata"] == {}
