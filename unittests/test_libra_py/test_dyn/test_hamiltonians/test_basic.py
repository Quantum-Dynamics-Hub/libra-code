import numpy as np

from libra_py.dyn.hamiltonians import (
    copy_hamiltonian_content,
    hamiltonian_memory_status,
    reset_hamiltonian_storage,
)


def test_hamiltonian_memory_status_and_reset(hamiltonian_storage):
    storage = hamiltonian_storage
    storage.ham_adi[0, 0] = np.diag([0.1, 0.2])
    storage.nac_adi[0, 0, 0, 1] = 0.03
    storage.ordering_adi[0, 0] = [1, 0]

    status = hamiltonian_memory_status(storage, der_lvl=2)

    assert status["ham_adi"] == 1
    assert status["d2ham_adi"] == 1

    reset_hamiltonian_storage(storage, der_lvl=2)

    assert np.allclose(storage.ham_adi[0, 0], 0.0)
    assert np.allclose(storage.nac_adi[0, 0], 0.0)
    assert np.allclose(storage.ordering_adi[0, 0], [0, 1])


def test_copy_hamiltonian_content(hamiltonian_storage):
    src = hamiltonian_storage
    dst = type(src)(
        backend=src.backend,
        ntraj=src.ntraj,
        ndof=src.ndof,
        nstates=src.nstates,
        ntbf_initial=src.ntbf_initial,
        ntbf_capacity=src.ntbf_capacity,
    )
    dst.allocate_hamiltonian_derivatives(der_lvl=2)

    src.ham_dia[0, 0] = np.array([[0.0, 0.1], [0.1, 0.2]])
    src.d1ham_dia[0, 0, 0] = np.diag([0.3, 0.4])
    src.gs_kinetic_energy = 2.5

    copy_hamiltonian_content(dst, src, der_lvl=2)

    assert np.allclose(dst.ham_dia, src.ham_dia)
    assert np.allclose(dst.d1ham_dia, src.d1ham_dia)
    assert dst.gs_kinetic_energy == 2.5
