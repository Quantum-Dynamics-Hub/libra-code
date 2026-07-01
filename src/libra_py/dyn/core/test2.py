"""
Demonstrate TensorStorage usage.

Run from this directory:

    python test2.py
"""

import numpy as np

import storage as st


def show(name, value):
    print(f"{name:32s}: {value}")


def main():
    storage = st.TensorStorage(
        backend=np,
        ntraj=3,
        ndof=2,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=4,
    )

    show("dimensions", (storage.ntraj, storage.ntbf, storage.ndof, storage.nstates))
    show("ndia/nadi/nnucl", (storage.ndia, storage.nadi, storage.nnucl))
    show("alive", storage.alive.astype(int).tolist())
    show("q shape", storage.q.shape)
    show("ampl_adi shape", storage.ampl_adi.shape)
    show("ham_adi shape", storage.ham_adi.shape)

    # Core nuclear variables.
    traj = 1
    tbf = 0
    storage.q[traj, tbf] = np.array([0.25, -0.10])
    storage.p[traj, tbf] = np.array([1.50, 0.00])
    storage.iM[traj, tbf] = np.array([1.0, 0.5])

    # Core electronic variables: adiabatic and diabatic amplitudes coexist.
    storage.ampl_adi[traj, tbf, 0] = 1.0
    storage.ampl_dia[traj, tbf, 1] = 1.0
    storage.dm_adi[traj, tbf] = np.outer(
        storage.ampl_adi[traj, tbf],
        storage.ampl_adi[traj, tbf].conjugate(),
    )

    # Core Hamiltonian variables.
    storage.ham_adi[traj, tbf] = np.diag([0.0, 0.1, 0.2])
    storage.nac_adi[traj, tbf, 0, 1] = 0.01
    storage.hvib_adi[traj, tbf] = storage.ham_adi[traj, tbf] - 1j * storage.nac_adi[traj, tbf]
    storage.act_states[traj, tbf] = 0

    show("q[1,0]", storage.q[traj, tbf].tolist())
    show("ampl_adi[1,0]", storage.ampl_adi[traj, tbf].tolist())
    show("ampl_dia[1,0]", storage.ampl_dia[traj, tbf].tolist())
    show("active state", int(storage.act_states[traj, tbf]))

    # Spawning adds a new TBF slot and marks it alive for one trajectory.
    new_tbf = storage.spawn(traj)
    storage.q[traj, new_tbf] = storage.q[traj, tbf]
    storage.p[traj, new_tbf] = storage.p[traj, tbf]
    storage.ampl_adi[traj, new_tbf, 2] = 1.0

    show("spawned tbf", new_tbf)
    show("ntbf after spawn", storage.ntbf)
    show("alive after spawn", storage.alive.astype(int).tolist())

    # Method-specific arrays are allocated only when requested.
    show("fssh2 status before", storage.fssh2_vars_status)
    show("dm_adi_prev before", storage.dm_adi_prev)
    storage.allocate_fssh2()
    storage.dm_adi_prev[traj, tbf] = storage.dm_adi[traj, tbf]
    show("fssh2 status after", storage.fssh2_vars_status)
    show("dm_adi_prev shape", storage.dm_adi_prev.shape)

    show("afssh status before", storage.afssh_vars_status)
    storage.allocate_afssh()
    storage.dR[traj, tbf, 0, 0, 0] = 0.001
    show("afssh status after", storage.afssh_vars_status)
    show("dR shape", storage.dR.shape)

    show("dc1_adi before", storage.dc1_adi)
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    storage.dc1_adi[traj, tbf, 0, 0, 1] = 0.02
    show("dc1_adi status", storage.dc1_adi_mem_status)
    show("dc1_adi shape", storage.dc1_adi.shape)
    show("d2ham_adi status", storage.d2ham_adi_mem_status)

    storage.allocate_simple_decoherence()
    storage.coherence_clocks[traj, tbf, 0, 1] = 5.0
    show("simple decoherence status", storage.simple_decoherence_vars_status)
    show("coherence_clocks shape", storage.coherence_clocks.shape)

    # Prune and compact remove slots that are inactive for all trajectories.
    storage.prune(traj, new_tbf)
    storage.compact()
    show("ntbf after prune/compact", storage.ntbf)
    show("q shape after compact", storage.q.shape)
    show("dR shape after compact", storage.dR.shape)

    assert storage.electronic_vars_status == 1
    assert storage.nuclear_vars_status == 1
    assert storage.fssh2_vars_status == 1
    assert storage.afssh_vars_status == 1
    assert storage.simple_decoherence_vars_status == 1
    assert storage.ampl_adi.shape[1] == storage.ntbf_capacity

    print("TensorStorage demonstration completed successfully.")


if __name__ == "__main__":
    main()
