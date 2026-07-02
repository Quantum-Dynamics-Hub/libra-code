import numpy as np

from libra_py.dyn.hamiltonians import HamiltonianEngine


def test_hamiltonian_engine_evaluate_writes_storage_and_returns_storage(
    hamiltonian_storage,
    trajectory,
    adiabatic_model_fn,
):
    storage = hamiltonian_storage
    engine = HamiltonianEngine(storage.backend)

    result = engine.evaluate(
        trajectory,
        storage,
        adiabatic_model_fn,
        rep="adiabatic",
    )

    assert result is storage
    assert np.allclose(storage.ham_adi[0, 0], np.diag([0.0, 0.2]))
    assert np.allclose(
        engine.active_matrix(storage, trajectory, rep="adiabatic", kind="vibronic"),
        storage.hvib_adi[0, [0]],
    )


def test_hamiltonian_engine_apply_ld_rotates_hamiltonian_and_nac(
    hamiltonian_storage,
    trajectory,
):
    storage = hamiltonian_storage
    engine = HamiltonianEngine(storage.backend)
    storage.ham_adi[0, 0] = np.diag([0.0, 0.2])
    storage.nac_adi[0, 0] = np.array([[0.0, 0.03], [-0.03, 0.0]])
    theta = np.pi / 4.0
    transform = np.array(
        [
            [
                [np.cos(theta), -np.sin(theta)],
                [np.sin(theta), np.cos(theta)],
            ]
        ],
        dtype=complex,
    )

    expected_ham = np.einsum(
        "...ij,...jk,...kl->...il",
        np.swapaxes(np.conjugate(transform), -1, -2),
        storage.ham_adi[0, [0]],
        transform,
    )

    engine.apply_ld(storage, trajectory, transform, reps=("adiabatic",))
    engine.build_hvib(storage, trajectory, reps=("adiabatic",))

    assert np.allclose(storage.ham_adi[0, [0]], expected_ham)
    assert np.allclose(storage.basis_transform[0, [0]], transform)
    assert np.allclose(
        storage.hvib_adi[0, 0],
        storage.ham_adi[0, 0] - 1j * storage.nac_adi[0, 0],
    )
