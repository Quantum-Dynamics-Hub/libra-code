"""
Educational example: TensorStorage.

TensorStorage owns the numerical arrays used by the Python dynamics prototype.
The first two axes are always:

    (ntraj, ntbf_capacity, ...)

where:

    ntraj          number of independent trajectories
    ntbf_capacity allocated number of TBF storage slots
    ndof           nuclear degrees of freedom
    nstates        electronic states

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/core/example_tensor_storage.py
"""

import numpy as np

from libra_py.dyn.core.storage import TensorStorage


def section(title):
    print(f"\n--- {title} ---")


def main():
    section("Create storage")

    # backend can be numpy-like. The storage class calls backend.zeros(...)
    # and keeps the door open for torch/jax-like backends later.
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=3,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=4,
    )

    print("ndia, nadi, nnucl:", storage.ndia, storage.nadi, storage.nnucl)
    print("q shape:", storage.q.shape)
    print("ampl_adi shape:", storage.ampl_adi.shape)
    print("ham_adi shape:", storage.ham_adi.shape)
    print("alive mask:\n", storage.alive.astype(int))

    section("Write core fields")

    traj = 0
    tbf = 0

    # Nuclear variables use C++ names from dyn_variables.h.
    storage.q[traj, tbf] = [0.0, 0.1, 0.2]
    storage.p[traj, tbf] = [1.0, 0.0, -1.0]
    storage.iM[traj, tbf] = [1.0, 0.5, 0.25]

    # Adiabatic and diabatic amplitudes coexist.
    storage.ampl_adi[traj, tbf] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.ampl_dia[traj, tbf] = [0.4 + 0.0j, 0.6 + 0.0j]

    # Hamiltonian-like data use nHamiltonian-style names.
    storage.ham_adi[traj, tbf] = np.diag([0.0, 0.25])
    storage.nac_adi[traj, tbf, 0, 1] = 0.01
    storage.hvib_adi[traj, tbf] = (
        storage.ham_adi[traj, tbf] - 1j * storage.nac_adi[traj, tbf]
    )

    print("q:", storage.q[traj, tbf])
    print("ampl_adi:", storage.ampl_adi[traj, tbf])
    print("ampl_dia:", storage.ampl_dia[traj, tbf])
    print("hvib_adi:\n", storage.hvib_adi[traj, tbf])

    section("Spawn, prune, compact")

    # spawn() activates the next global TBF slot for one trajectory.
    spawned_slot = storage.spawn(traj)
    storage.q[traj, spawned_slot] = storage.q[traj, tbf]
    storage.ampl_adi[traj, spawned_slot, 1] = 1.0

    print("spawned slot:", spawned_slot)
    print("ntbf:", storage.ntbf)
    print("alive mask:\n", storage.alive.astype(int))

    # prune() marks a trajectory/TBF slot inactive. compact() removes slots
    # that are inactive for all trajectories.
    storage.prune(traj, spawned_slot)
    storage.compact()
    print("q shape after compact:", storage.q.shape)

    section("Method-specific allocation")

    # Large or method-specific arrays are not allocated by default.
    print("dm_adi_prev before allocate_fssh2:", storage.dm_adi_prev)
    storage.allocate_fssh2()
    storage.dm_adi_prev[traj, tbf] = np.eye(storage.nstates)
    print("fssh2 status:", storage.fssh2_vars_status)
    print("dm_adi_prev shape:", storage.dm_adi_prev.shape)


if __name__ == "__main__":
    main()
