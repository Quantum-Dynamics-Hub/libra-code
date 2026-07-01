"""
Educational example: TrajectoryBasisFunction.

A TrajectoryBasisFunction object stores topology and metadata for one nuclear
wavepacket. The numerical state lives in TensorStorage and is reached through
the TBF's TensorView.

Important ids:

    tbf.id     object/registry id
    tbf.tbf_id TensorStorage slot id

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/core/example_tbf.py
"""

import numpy as np

from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.tbf import TrajectoryBasisFunction
from libra_py.dyn.core.views import TensorView


def section(title):
    print(f"\n--- {title} ---")


def main():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=2,
        ntbf_initial=1,
    )

    view = TensorView(storage, traj_id=1, tbf_id=0)
    tbf = TrajectoryBasisFunction(
        id=10,
        view=view,
        trajectory_id=1,
        metadata={"label": "initial"},
    )

    section("Identity and topology")

    print(tbf)
    print("object id:", tbf.id)
    print("storage slot id:", tbf.tbf_id)
    print("trajectory id:", tbf.trajectory_id)

    section("Write data through the view")

    tbf.view.q[:] = [0.0, 0.5]
    tbf.view.p[:] = [1.0, -1.0]
    tbf.view.ampl_adi[:] = [1.0 + 0.0j, 0.0 + 0.0j]

    print("q through view:", tbf.view.q)
    print("same q in storage:", storage.q[1, 0])

    section("Lifecycle")

    tbf.deactivate()
    print("object alive:", tbf.alive)
    print("storage alive:", storage.alive[1, 0])

    tbf.activate()
    tbf.set_parent(parent_id=3)
    tbf.set_spawn_time(time=4.0)

    print("summary:", tbf.to_dict())


if __name__ == "__main__":
    main()
