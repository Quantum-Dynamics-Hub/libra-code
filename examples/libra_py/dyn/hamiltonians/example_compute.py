"""
Educational example: Hamiltonian construction.

This example evaluates a tiny two-state diabatic model, transforms it to the
adiabatic representation, and shows where the tensors live in TensorStorage.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/hamiltonians/example_compute.py
"""

import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import compute_adiabatic, compute_diabatic


def section(title):
    print(f"\n--- {title} ---")


def model(R, P, storage, traj, params=None):
    return {
        "ham_dia": np.array([[[0.0, 0.01], [0.01, 0.1]]], dtype=complex),
        "ovlp_dia": np.array([np.eye(2)], dtype=complex),
        "nac_dia": np.zeros((1, 2, 2), dtype=complex),
        "d1ham_dia": np.array([[[[0.1, 0.0], [0.0, 0.2]]]], dtype=complex),
        "dc1_dia": np.zeros((1, 1, 2, 2), dtype=complex),
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
    storage.q[0, 0] = [0.25]
    storage.p[0, 0] = [0.5]

    section("Diabatic model")
    compute_diabatic(storage, traj, model, der_lvl=1)
    print("ham_dia:\n", storage.ham_dia[0, 0])
    print("hvib_dia:\n", storage.hvib_dia[0, 0])

    section("Adiabatic transform")
    compute_adiabatic(storage, traj, from_diabatic=True, der_lvl=1)
    print("ham_adi:\n", storage.ham_adi[0, 0])
    print("basis_transform:\n", storage.basis_transform[0, 0])
    print("d1ham_adi:\n", storage.d1ham_adi[0, 0])


if __name__ == "__main__":
    main()
