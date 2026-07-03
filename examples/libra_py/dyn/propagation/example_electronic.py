"""
Educational example: TD-SE electronic propagation.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/propagation/example_electronic.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.propagation import tdse_step, update_density


def section(title):
    print(f"\n--- {title} ---")


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

    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.hvib_adi[0, 0] = np.array(
        [[0.0, 0.02j], [-0.02j, 0.2]],
        dtype=complex,
    )

    section("Initial amplitudes")
    print("C:", storage.ampl_adi[0, 0])

    section("TD-SE step")
    tdse_step(traj, storage, dt=0.5, rep="adiabatic", hamiltonian_type="vibronic")
    update_density(storage, traj, rep="adiabatic")
    print("C:", storage.ampl_adi[0, 0])
    print("rho:\n", storage.dm_adi[0, 0])


if __name__ == "__main__":
    main()
