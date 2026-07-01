"""
Demonstrate TrajectoryBasisFunction usage.

Run from this directory:

    python test4.py
"""

from pathlib import Path
import sys

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[3]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.tbf import TrajectoryBasisFunction
from libra_py.dyn.core.views import TensorView


def show(name, value):
    print(f"{name:34s}: {value}")


def main():
    # TensorStorage owns numerical arrays; TensorView points to one
    # (trajectory, TBF-slot) location; TrajectoryBasisFunction owns topology
    # metadata and lifecycle flags.
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=3,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=4,
    )

    view = TensorView(storage, traj_id=1, tbf_id=0)
    tbf = TrajectoryBasisFunction(
        id=7,
        view=view,
        trajectory_id=1,
        parent_id=None,
        spawn_time=0.0,
        metadata={"label": "initial packet"},
    )

    show("tbf repr", tbf)
    show("object id", tbf.object_id)
    show("storage tbf slot", tbf.tbf_id)
    show("trajectory id", tbf.trajectory_id)
    show("alive in storage", tbf.is_alive_in_storage())

    # The TBF view is a writable window into TensorStorage.
    tbf.view.q[:] = [0.0, 0.1, 0.2]
    tbf.view.p[:] = [1.0, 0.0, -1.0]
    tbf.view.ampl_adi[:] = [1.0 + 0.0j, 0.0 + 0.0j]
    tbf.view.ham_adi[:] = np.diag([0.0, 0.5])

    show("view q", tbf.view.q.tolist())
    show("storage q[1,0]", storage.q[1, 0].tolist())
    show("view ampl_adi", tbf.view.ampl_adi.tolist())

    # Lifecycle methods update both the object flag and TensorStorage.alive by
    # default. This keeps object bookkeeping and tensor bookkeeping aligned.
    tbf.deactivate()
    show("alive after deactivate", tbf.alive)
    show("storage alive[1,0]", bool(storage.alive[1, 0]))

    tbf.activate()
    show("alive after activate", tbf.alive)
    show("storage alive[1,0]", bool(storage.alive[1, 0]))

    # Parent and spawn-time metadata can be adjusted independently of storage.
    tbf.set_parent(parent_id=3)
    tbf.set_spawn_time(time=12.5)
    tbf.metadata["reason"] = "population threshold"

    show("parent id", tbf.parent_id)
    show("spawn time", tbf.spawn_time)
    show("metadata", tbf.metadata)
    show("summary", tbf.to_dict())

    assert tbf.id == 7
    assert tbf.tbf_id == 0
    assert tbf.traj_id == 1
    assert storage.q[1, 0, 2] == 0.2
    assert tbf.is_alive_in_storage()

    print("TrajectoryBasisFunction demonstration completed successfully.")


if __name__ == "__main__":
    main()
