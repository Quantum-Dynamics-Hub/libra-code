"""
Educational example: Trajectory.

Trajectory represents one electron-nuclear wavefunction. It stores TBF storage
slot ids in trajectory.tbf_ids, so the list can be used directly for tensor
indexing:

    storage.q[trajectory.id, trajectory.tbf_ids]

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/core/example_trajectory.py
"""

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


def section(title):
    print(f"\n--- {title} ---")


def main():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=4,
    )
    traj = Trajectory(id=1)
    traj.metadata["label"] = "demo trajectory"

    section("Attach TBF objects")

    tbf0 = make_tbf(storage, object_id=20, traj_id=1, slot=0)
    traj.add_tbf(tbf0)
    tbf0.view.q[:] = [0.0, 0.1]

    slot1 = storage.spawn(traj.id)
    tbf1 = make_tbf(storage, object_id=21, traj_id=1, slot=slot1)
    traj.add_tbf(tbf1)
    tbf1.view.q[:] = [0.5, 0.6]

    print("trajectory:", traj)
    print("storage slot ids:", traj.tbf_ids)
    print("object ids:", traj.tbf_object_ids)
    print("slot -> object:", traj.tbf_object_by_slot)

    section("Use tbf_ids for tensor slicing")

    print("q for trajectory:\n", storage.q[traj.id, traj.tbf_ids])

    section("Lifecycle bookkeeping")

    tbf1.deactivate()
    print("active slots:", traj.active_tbf_ids(storage))

    traj.remove_tbf(slot1)
    traj.set_norm(0.99)
    print("summary:", traj.to_dict())


if __name__ == "__main__":
    main()
