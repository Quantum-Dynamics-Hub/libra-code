"""
Educational example: compact dynamics observables.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/observables/example_observables.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.observables import ObservableConfig, compute_observables


def main():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
    )
    traj = Trajectory(0)
    traj.tbf_ids = [0]

    storage.p[0, 0] = [1.0]
    storage.iM[0, 0] = [0.5]
    storage.ampl_adi[0, 0] = [0.6 + 0.0j, 0.8 + 0.0j]
    storage.act_states[0, 0] = 1
    storage.ham_adi[0, 0] = np.diag([0.0, 0.2])

    data = compute_observables(
        storage,
        traj,
        ObservableConfig(include_coordinates=True),
        step=storage.timestep,
        time=0.0,
    )

    for key, value in data.items():
        print(key, "=", value)

    print("\nold save.py-style names")
    legacy = compute_observables(
        storage,
        traj,
        ObservableConfig(output_level=2),
        step=storage.timestep,
        time=0.0,
    )
    print("Ekin_ave =", legacy["Ekin_ave"])
    print("se_pop_adi =", legacy["se_pop_adi"])
    print("sh_pop_adi =", legacy["sh_pop_adi"])


if __name__ == "__main__":
    main()
