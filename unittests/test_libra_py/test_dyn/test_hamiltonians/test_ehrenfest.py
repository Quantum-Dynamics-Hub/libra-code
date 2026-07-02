import numpy as np

from libra_py.dyn.hamiltonians import (
    compute_adiabatic,
    compute_diabatic,
    ehrenfest_energy_adi,
    ehrenfest_energy_dia,
    ehrenfest_force_tensors_adi,
    ehrenfest_force_tensors_dia,
    ehrenfest_forces_adi,
    ehrenfest_forces_dia,
)


def test_ehrenfest_energies(
    hamiltonian_storage,
    trajectory,
    diabatic_model_fn,
    adiabatic_model_fn,
):
    storage = hamiltonian_storage
    compute_diabatic(storage, trajectory, diabatic_model_fn, der_lvl=1)
    compute_adiabatic(storage, trajectory, adiabatic_model_fn, der_lvl=1)
    storage.ampl_dia[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.ampl_adi[0, 0] = [0.0 + 0.0j, 1.0 + 0.0j]

    assert np.allclose(ehrenfest_energy_dia(storage, trajectory), [0.0])
    assert np.allclose(ehrenfest_energy_adi(storage, trajectory), [0.2])


def test_ehrenfest_force_tensors_and_forces(
    hamiltonian_storage,
    trajectory,
    diabatic_model_fn,
    adiabatic_model_fn,
):
    storage = hamiltonian_storage
    compute_diabatic(storage, trajectory, diabatic_model_fn, der_lvl=1)
    compute_adiabatic(storage, trajectory, adiabatic_model_fn, der_lvl=1)
    storage.ampl_dia[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]
    storage.ampl_adi[0, 0] = [1.0 + 0.0j, 0.0 + 0.0j]

    tensors_dia = ehrenfest_force_tensors_dia(storage, trajectory, option=1)
    tensors_adi = ehrenfest_force_tensors_adi(storage, trajectory, option=1)
    forces_dia = ehrenfest_forces_dia(storage, trajectory, option=1)
    forces_adi = ehrenfest_forces_adi(storage, trajectory, option=1)

    assert tensors_dia.shape == (1, 1, 2, 2)
    assert tensors_adi.shape == (1, 1, 2, 2)
    assert np.allclose(forces_dia, [[-0.1]])
    assert np.allclose(forces_adi, [[-0.1]])


def test_ehrenfest_adi_gamma_term_is_normalized(
    hamiltonian_storage,
    trajectory,
    adiabatic_model_fn,
):
    storage = hamiltonian_storage
    compute_adiabatic(storage, trajectory, adiabatic_model_fn, der_lvl=1)
    storage.ampl_adi[0, 0] = [2.0 + 0.0j, 0.0 + 0.0j]

    forces = ehrenfest_forces_adi(
        storage,
        trajectory,
        option=1,
        gamma=0.5,
    )

    # Base force: -<C|dH|C>/norm = -0.1.
    # Gamma term: gamma * trace(dH)/norm = 0.5 * 0.3 / 4.
    assert np.allclose(forces, [[-0.0625]])
