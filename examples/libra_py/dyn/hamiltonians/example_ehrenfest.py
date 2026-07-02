"""
Educational example: Ehrenfest energies and forces.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/hamiltonians/example_ehrenfest.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import (
    compute_adiabatic,
    compute_diabatic,
    ehrenfest_energy_adi,
    ehrenfest_energy_dia,
    ehrenfest_forces_adi,
    ehrenfest_forces_dia,
)


def section(title):
    print(f"\n--- {title} ---")


def diabatic_model(R, P, storage, traj, params=None):
    return {
        "ham_dia": np.array([[[0.0, 0.01], [0.01, 0.1]]], dtype=complex),
        "ovlp_dia": np.array([np.eye(2)], dtype=complex),
        "d1ham_dia": np.array([[[[0.1, 0.0], [0.0, 0.2]]]], dtype=complex),
        "dc1_dia": np.zeros((1, 1, 2, 2), dtype=complex),
    }


def adiabatic_model(R, P, storage, traj, params=None):
    return {
        "H_adi": np.array([np.diag([0.0, 0.2])], dtype=complex),
        "dH_adi": np.array([[[[0.1, 0.0], [0.0, 0.2]]]], dtype=complex),
        "DC1_adi": np.zeros((1, 1, 2, 2), dtype=complex),
    }


def main():
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=1,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)

    traj = Trajectory(0)
    traj.tbf_ids = [0]
    storage.ampl_dia[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.ampl_adi[0, 0] = [0.0 + 0.0j, 1.0 + 0.0j]

    compute_diabatic(storage, traj, diabatic_model, der_lvl=1)
    compute_adiabatic(storage, traj, adiabatic_model, der_lvl=1)

    section("Energies")
    print("E_dia:", ehrenfest_energy_dia(storage, traj))
    print("E_adi:", ehrenfest_energy_adi(storage, traj))

    section("Gradient-only forces")
    print("F_dia:", ehrenfest_forces_dia(storage, traj, option=1))
    print("F_adi:", ehrenfest_forces_adi(storage, traj, option=1))


if __name__ == "__main__":
    main()
