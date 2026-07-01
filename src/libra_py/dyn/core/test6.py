"""
Demonstrate TrajectoryManager usage.

Run from this directory:

    python test6.py
"""

from pathlib import Path
import sys

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[3]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from libra_py.dyn.core.manager import TrajectoryManager
from libra_py.dyn.core.storage import TensorStorage


def show(name, value):
    print(f"{name:36s}: {value}")


def main():
    # The manager owns topology objects and the shared TensorStorage reference.
    storage = TensorStorage(
        backend=np,
        ntraj=3,
        ndof=2,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=5,
    )
    manager = TrajectoryManager(storage)

    # Create two independent trajectory objects. Their ids match the first
    # TensorStorage index used by their TBF views.
    traj0 = manager.create_trajectory(metadata={"label": "left branch"})
    traj1 = manager.create_trajectory(metadata={"label": "right branch"})

    show("trajectory ids", (traj0.id, traj1.id))

    # Register initial TBF objects for storage slot 0 on each trajectory.
    tbf0 = manager.create_tbf(traj0, idx=0, metadata={"role": "initial"})
    tbf1 = manager.create_tbf(traj1, idx=0, metadata={"role": "initial"})

    # Fill data through TBF views. The data live in storage, not in the TBF
    # objects themselves.
    tbf0.view.q[:] = [0.0, 0.5]
    tbf0.view.p[:] = [1.0, 0.0]
    tbf0.view.ampl_adi[:] = [1.0, 0.0, 0.0]

    tbf1.view.q[:] = [2.0, 2.5]
    tbf1.view.p[:] = [0.0, -1.0]
    tbf1.view.ampl_adi[:] = [0.0, 1.0, 0.0]

    show("manager summary counts", {
        "n_trajectories": manager.summary()["n_trajectories"],
        "n_tbfs": manager.summary()["n_tbfs"],
    })
    show("traj0 slots", traj0.tbf_ids)
    show("traj1 slots", traj1.tbf_ids)

    # Spawn a new TBF on traj0. With parent=tbf0 and copy_from_parent=True,
    # manager.copy_tbf_data() copies core fields into the new storage slot.
    spawned = manager.spawn_tbf(
        traj0,
        parent=tbf0,
        spawn_time=10.0,
        metadata={"role": "spawned"},
    )

    show("spawned object id", spawned.id)
    show("spawned storage slot", spawned.tbf_id)
    show("spawned parent id", spawned.parent_id)
    show("traj0 slots after spawn", traj0.tbf_ids)
    show("spawned q copied", spawned.view.q.tolist())
    show("spawned ampl_adi copied", spawned.view.ampl_adi.tolist())

    # Modify the spawned packet independently after copying.
    spawned.view.q[:] += [0.1, -0.1]
    spawned.view.ampl_adi[:] = [0.0, 0.0, 1.0]

    show("parent q unchanged", tbf0.view.q.tolist())
    show("spawned q modified", spawned.view.q.tolist())
    show("active traj0 slots", manager.active_tbf_ids(traj0))
    show("active traj0 objects", [tbf.id for tbf in manager.active_tbfs(traj0)])

    # Deactivate/reactivate a TBF. This updates both the object and
    # TensorStorage.alive bookkeeping.
    manager.deactivate_tbf(spawned)
    show("active after deactivate", manager.active_tbf_ids(traj0))
    show("storage alive spawned", bool(storage.alive[traj0.id, spawned.tbf_id]))

    manager.activate_tbf(spawned.id)
    show("active after reactivate", manager.active_tbf_ids(traj0))

    # Deactivate a whole trajectory and its TBFs.
    manager.deactivate_trajectory(traj1, deactivate_tbfs=True)
    show("traj1 active", traj1.active)
    show("traj1 active slots", manager.active_tbf_ids(traj1))

    manager.activate_trajectory(traj1, activate_tbfs=True)
    show("traj1 active after activate", traj1.active)
    show("traj1 active slots restored", manager.active_tbf_ids(traj1))

    # Remove the spawned TBF object. The storage slot is deactivated and the
    # object is removed from the manager registry.
    removed = manager.remove_tbf(spawned.id)
    show("removed object id", removed.id)
    show("spawned in registry", spawned.id in manager.tbfs)
    show("traj0 slots after remove", traj0.tbf_ids)

    # compact_storage() drops globally inactive storage slots and remaps any
    # surviving TBF views. Object ids are preserved; storage slot ids may change.
    slot_map = manager.compact_storage()
    show("compact slot map", slot_map)
    show("storage ntbf_capacity", storage.ntbf_capacity)
    show("tbf0 slot after compact", tbf0.tbf_id)
    show("tbf1 slot after compact", tbf1.tbf_id)

    # The manager summary is a compact serializable snapshot of topology.
    summary = manager.summary()
    show("summary n_tbfs", summary["n_tbfs"])
    show("summary traj0", summary["trajectories"][traj0.id])
    show("summary tbf0", summary["tbfs"][tbf0.id])

    assert summary["n_trajectories"] == 2
    assert summary["n_tbfs"] == 2
    assert tbf0.id in manager.tbfs
    assert tbf1.id in manager.tbfs
    assert spawned.id not in manager.tbfs
    assert manager.active_tbf_ids(traj0) == [tbf0.tbf_id]
    assert manager.active_tbf_ids(traj1) == [tbf1.tbf_id]

    print("TrajectoryManager demonstration completed successfully.")


if __name__ == "__main__":
    main()
