"""
Educational example: coupled propagation helpers.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/propagation/example_coupled.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.propagation import (
    ehrenfest_forces,
    state_specific_forces,
    update_density,
)


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
    storage.allocate_hamiltonian_derivatives()

    traj = Trajectory(0)
    traj.tbf_ids = [0]

    storage.act_states[0, 0] = 1
    storage.ampl_adi[0, 0] = [0.6 + 0.0j, 0.8 + 0.0j]
    storage.ham_adi[0, 0] = np.diag([0.0, 0.2])
    storage.d1ham_adi[0, 0] = np.array(
        [
            [[0.1, 0.0], [0.0, 0.3]],
            [[-0.2, 0.0], [0.0, 0.4]],
        ],
        dtype=complex,
    )
    storage.dc1_adi[0, 0] = np.zeros((2, 2, 2), dtype=complex)

    section("Density")
    print("rho:\n", update_density(storage, traj, rep="adiabatic")[0])

    section("Forces")
    print("state-specific:", state_specific_forces(storage, traj))
    print("Ehrenfest:", ehrenfest_forces(storage, traj, option=1))


if __name__ == "__main__":
    main()
