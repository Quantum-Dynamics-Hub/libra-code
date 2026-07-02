import numpy as np

from libra_py.dyn.hamiltonians import (
    compute_adiabatic,
    compute_adiabatic_from_diabatic,
    compute_diabatic,
)
from libra_py.dyn.utils.linalg import generalized_eigh


def test_compute_diabatic_writes_model_fields_and_builds_hvib(
    hamiltonian_storage,
    trajectory,
    diabatic_model_fn,
):
    storage = hamiltonian_storage

    compute_diabatic(storage, trajectory, diabatic_model_fn, der_lvl=2)

    assert np.allclose(storage.ham_dia[0, 0], [[0.0, 0.01], [0.01, 0.1]])
    assert np.allclose(storage.ovlp_dia[0, 0], np.eye(2))
    assert np.allclose(
        storage.hvib_dia[0, 0],
        storage.ham_dia[0, 0] - 1j * storage.nac_dia[0, 0],
    )
    assert storage.d2ham_dia.shape == (1, 2, 1, 1, 2, 2)


def test_compute_adiabatic_from_model_aliases_and_builds_hvib(
    hamiltonian_storage,
    trajectory,
    adiabatic_model_fn,
):
    storage = hamiltonian_storage

    compute_adiabatic(storage, trajectory, adiabatic_model_fn, der_lvl=1)

    assert np.allclose(storage.ham_adi[0, 0], np.diag([0.0, 0.2]))
    assert np.allclose(
        storage.hvib_adi[0, 0],
        storage.ham_adi[0, 0] - 1j * storage.nac_adi[0, 0],
    )
    assert np.allclose(storage.d1ham_adi[0, 0, 0], np.diag([0.1, 0.2]))


def test_compute_adiabatic_from_diabatic_diagonalizes_and_transforms_derivatives(
    hamiltonian_storage,
    trajectory,
    diabatic_model_fn,
):
    storage = hamiltonian_storage
    compute_diabatic(storage, trajectory, diabatic_model_fn, der_lvl=1)

    compute_adiabatic_from_diabatic(storage, trajectory, der_lvl=1)

    expected_energy, expected_transform = generalized_eigh(
        storage.ham_dia[0, [0]],
        storage.ovlp_dia[0, [0]],
    )
    expected_ham = np.zeros((2, 2), dtype=complex)
    expected_ham[np.diag_indices(2)] = expected_energy[0]

    assert np.allclose(storage.ham_adi[0, 0], expected_ham)
    assert np.allclose(storage.basis_transform[0, 0], expected_transform[0])
    assert np.allclose(storage.hvib_adi[0, 0], storage.ham_adi[0, 0])
    assert np.allclose(storage.d1ham_adi[0, 0, 0], np.diag(np.diag(storage.d1ham_adi[0, 0, 0])))
    assert np.allclose(storage.dc1_adi[0, 0, 0], -storage.dc1_adi[0, 0, 0].T)
