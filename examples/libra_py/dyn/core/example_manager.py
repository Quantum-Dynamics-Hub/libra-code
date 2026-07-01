"""
Educational example: TrajectoryManager.

TrajectoryManager owns the topology registries and coordinates TensorStorage,
Trajectory, TrajectoryBasisFunction, and TensorView.

Use it when object ids and storage slot ids must stay consistent.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/core/example_manager.py
"""

import numpy as np

from libra_py.dyn.core.manager import TrajectoryManager
from libra_py.dyn.core.storage import TensorStorage


def section(title):
    print(f"\n--- {title} ---")


def main():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=4,
    )
    manager = TrajectoryManager(storage)

    section("Create topology")

    traj = manager.create_trajectory(metadata={"label": "demo"})
    parent = manager.create_tbf(traj, idx=0, metadata={"role": "parent"})

    parent.view.q[:] = [0.0, 0.5]
    parent.view.p[:] = [1.0, 0.0]
    parent.view.ampl_adi[:] = [1.0, 0.0, 0.0]

    print("trajectory id:", traj.id)
    print("parent object id:", parent.id)
    print("parent storage slot:", parent.tbf_id)

    section("Spawn from a parent")

    child = manager.spawn_tbf(
        traj,
        parent=parent,
        spawn_time=2.5,
        metadata={"role": "spawned"},
    )

    # The manager copied core fields from parent to child. We can then modify
    # the child independently.
    child.view.q[:] += [0.1, -0.1]
    child.view.ampl_adi[:] = [0.0, 0.0, 1.0]

    print("active slots:", manager.active_tbf_ids(traj))
    print("parent q:", parent.view.q)
    print("child q:", child.view.q)
    print("child parent id:", child.parent_id)

    section("Lifecycle")

    manager.deactivate_tbf(child)
    print("active after child deactivate:", manager.active_tbf_ids(traj))

    manager.activate_tbf(child)
    print("active after child activate:", manager.active_tbf_ids(traj))

    removed = manager.remove_tbf(child.id)
    print("removed object id:", removed.id)
    print("registered TBF ids:", sorted(manager.tbfs))

    section("Compact storage")

    # Compaction drops globally inactive storage slots and remaps surviving
    # TBF views. Object ids remain stable.
    slot_map = manager.compact_storage()
    print("slot map after compact:", slot_map)
    print("parent storage slot after compact:", parent.tbf_id)

    section("Summary")

    print(manager.summary())


if __name__ == "__main__":
    main()
