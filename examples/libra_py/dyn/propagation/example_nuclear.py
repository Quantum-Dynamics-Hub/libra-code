"""
Educational example: nuclear propagation kernels.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/propagation/example_nuclear.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.propagation import drift, kick, velocity


def section(title):
    print(f"\n--- {title} ---")


def main():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=2,
        nstates=2,
        ntbf_initial=1,
    )
    traj = Trajectory(0)
    traj.tbf_ids = [0]

    storage.q[0, 0] = [0.0, 1.0]
    storage.p[0, 0] = [2.0, -1.0]
    storage.iM[0, 0] = [0.5, 0.25]
    storage.f[0, 0] = [-0.2, 0.4]

    section("Velocity")
    print("v:", velocity(storage, traj))

    section("Drift")
    drift(storage, traj, dt=0.1)
    print("q:", storage.q[0, 0])

    section("Kick")
    kick(storage, traj, dt=0.1)
    print("p:", storage.p[0, 0])


if __name__ == "__main__":
    main()
