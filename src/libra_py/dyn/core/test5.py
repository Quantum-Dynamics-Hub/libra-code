"""
Demonstrate Trajectory usage.

Run from this directory:

    python test5.py
"""

from pathlib import Path
import sys

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[3]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.tbf import TrajectoryBasisFunction
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.core.views import TensorView


def show(name, value):
    print(f"{name:34s}: {value}")


def make_tbf(storage, object_id, traj_id, tbf_id, parent_id=None, spawn_time=None):
    view = TensorView(storage, traj_id=traj_id, tbf_id=tbf_id)
    return TrajectoryBasisFunction(
        id=object_id,
        view=view,
        trajectory_id=traj_id,
        parent_id=parent_id,
        spawn_time=spawn_time,
    )


def main():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=5,
    )

    traj = Trajectory(id=1)
    traj.metadata["ensemble"] = "demo"

    # Attach the initial TBF in storage slot 0.
    tbf0 = make_tbf(storage, object_id=10, traj_id=1, tbf_id=0)
    traj.add_tbf(tbf0)
    tbf0.view.q[:] = [0.1, 0.2]
    tbf0.view.ampl_adi[0] = 1.0

    show("trajectory after tbf0", traj)
    show("tbf storage slots", traj.tbf_ids)
    show("tbf object ids", traj.tbf_object_ids)

    # Spawn a new storage slot for this trajectory, then attach a second TBF
    # object to that slot. The object id and storage slot id are distinct.
    slot1 = storage.spawn(traj.id)
    tbf1 = make_tbf(
        storage,
        object_id=11,
        traj_id=1,
        tbf_id=slot1,
        parent_id=tbf0.id,
        spawn_time=5.0,
    )
    traj.add_tbf(tbf1)
    tbf1.view.q[:] = tbf0.view.q
    tbf1.view.ampl_adi[2] = 1.0

    show("spawned storage slot", slot1)
    show("trajectory after tbf1", traj)
    show("slot->object map", traj.tbf_object_by_slot)
    show("active slots", traj.active_tbf_ids(storage))

    # The trajectory stores storage slots in tbf_ids, so these values can be
    # used directly to slice TensorStorage arrays.
    q_active = storage.q[traj.id, traj.tbf_ids]
    ampl_active = storage.ampl_adi[traj.id, traj.tbf_ids]
    show("q for trajectory TBFs", q_active.tolist())
    show("ampl_adi for trajectory TBFs", ampl_active.tolist())

    # Deactivating a TBF updates TensorStorage.alive. The trajectory can then
    # report only storage-active TBF slots.
    tbf1.deactivate()
    show("active slots after tbf1 off", traj.active_tbf_ids(storage))
    show("all known slots", list(traj))

    # Removing a TBF removes the storage slot from trajectory bookkeeping and
    # also removes the corresponding object id from tbf_object_ids.
    traj.remove_tbf(slot1)
    show("trajectory after remove", traj)
    show("slot->object map after remove", traj.tbf_object_by_slot)

    # Wavefunction-level metadata lives on the trajectory, not on individual
    # storage tensors.
    traj.set_norm(0.998)
    traj.deactivate()
    show("trajectory active", traj.active)
    show("trajectory norm", traj.norm)
    show("summary", traj.to_dict())

    traj.activate()

    assert traj.has_tbf(0)
    assert not traj.has_tbf(slot1)
    assert len(traj) == 1
    assert traj.active
    assert traj.tbf_ids == [0]
    assert traj.tbf_object_ids == [10]
    assert traj.active_tbf_ids(storage) == [0]

    print("Trajectory demonstration completed successfully.")


if __name__ == "__main__":
    main()
