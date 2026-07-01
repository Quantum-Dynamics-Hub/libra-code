import numpy as np

from libra_py.dyn.core.storage import TensorStorage


def test_tensor_storage_core_allocation_and_fields():
    storage = TensorStorage(
        backend=np,
        ntraj=3,
        ndof=2,
        nstates=3,
        ntbf_initial=1,
        ntbf_capacity=4,
    )

    assert storage.ndia == 3
    assert storage.nadi == 3
    assert storage.nnucl == 2
    assert storage.q.shape == (3, 4, 2)
    assert storage.ampl_adi.shape == (3, 4, 3)
    assert storage.ham_adi.shape == (3, 4, 3, 3)
    assert storage.electronic_vars_status == 1
    assert storage.nuclear_vars_status == 1
    assert storage.afssh_vars_status == 0

    storage.q[1, 0] = [0.25, -0.1]
    storage.p[1, 0] = [1.5, 0.0]
    storage.ampl_adi[1, 0, 0] = 1.0
    storage.ampl_dia[1, 0, 1] = 1.0

    assert np.allclose(storage.q[1, 0], [0.25, -0.1])
    assert storage.ampl_adi[1, 0, 0] == 1.0
    assert storage.ampl_dia[1, 0, 1] == 1.0


def test_tensor_storage_spawn_prune_compact():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=2,
        nstates=2,
        ntbf_initial=1,
        ntbf_capacity=4,
    )

    slot = storage.spawn(traj=1)
    assert slot == 1
    assert storage.ntbf == 2
    assert bool(storage.alive[1, slot])

    storage.q[1, slot] = [1.0, 2.0]
    storage.prune(1, slot)
    assert not bool(storage.alive[1, slot])

    storage.compact()
    assert storage.ntbf == 1
    assert storage.ntbf_capacity == 1
    assert storage.q.shape == (2, 1, 2)


def test_tensor_storage_specific_allocators():
    storage = TensorStorage(
        backend=np,
        ntraj=2,
        ndof=3,
        nstates=2,
        ntbf_initial=1,
    )

    assert storage.dm_adi_prev is None
    assert storage.dR is None
    assert storage.dc1_adi is None
    assert storage.d2ham_adi is None

    storage.allocate_fssh2()
    storage.allocate_afssh()
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    storage.allocate_simple_decoherence()

    assert storage.fssh2_vars_status == 1
    assert storage.afssh_vars_status == 1
    assert storage.dc1_adi_mem_status == 1
    assert storage.d2ham_adi_mem_status == 0
    assert storage.dm_adi_prev.shape == (2, 16, 2, 2)
    assert storage.dR.shape == (2, 16, 3, 2, 2)
    assert storage.dc1_adi.shape == (2, 16, 3, 2, 2)
    assert storage.coherence_clocks.shape == (2, 16, 2, 2)

    storage.allocate_hamiltonian_derivatives(der_lvl=2)
    assert storage.d2ham_adi_mem_status == 1
    assert storage.d2ham_adi.shape == (2, 16, 3, 3, 2, 2)
