import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import HamiltonianEngine
from libra_py.dyn.models import (
    ShinMetiuDVRData,
    ShinMetiuDVRModel,
    bundled_shin_metiu_dvr_path,
    kinetic_energy_matrix,
    load_shin_metiu_dvr_data,
)


def test_kinetic_energy_matrix_matches_legacy_formula():
    matrix = kinetic_energy_matrix(5, 0.25)

    assert matrix.shape == (5, 5)
    assert np.allclose(matrix, matrix.T)
    assert np.all(matrix.diagonal() > 0.0)
    assert np.isclose(matrix[0, 0], 0.5 * np.pi**2 / (3.0 * 0.25**2) * (1.0 + 2.0 / 25.0))


def test_shin_metiu_dvr_data_interpolates_and_builds_covariant_derivative():
    data = _synthetic_dvr_data()
    result = data.interpolate(np.array([-0.5, 0.5]))

    expected_e0 = np.array([0.05, 0.15])
    expected_e1 = np.array([0.55, 0.65])
    expected_nac = np.array([0.005, 0.015])

    assert np.allclose(result["eigvals"][:, 0], expected_e0)
    assert np.allclose(result["eigvals"][:, 1], expected_e1)
    assert np.allclose(result["nac"][:, 0, 1], expected_nac)
    assert np.allclose(result["H_adi"][:, 0, 0], expected_e0)
    assert np.allclose(result["H_adi"][:, 1, 1], expected_e1)

    expected_dh01 = result["d_V"][:, 0, 1] + result["H_adi"][:, 0, 0] * expected_nac
    expected_dh01 -= expected_nac * result["H_adi"][:, 1, 1]
    assert np.allclose(result["dH_adi"][:, 0, 1], expected_dh01)


def test_shin_metiu_dvr_model_writes_adiabatic_quantities_through_engine():
    model = ShinMetiuDVRModel(dvr_data=_synthetic_dvr_data())
    q = np.array([-0.5, 0.5])

    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=2,
        ntbf_capacity=2,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    storage.q[0, :, 0] = q
    storage.p[0, :, 0] = [2.0, 3.0]
    storage.iM[0, :, 0] = 0.5

    traj = Trajectory(0)
    traj.tbf_ids = [0, 1]
    HamiltonianEngine(backend).evaluate(traj, storage, model, rep="adiabatic")

    direct = model.evaluate(np.asarray([q]))
    assert np.allclose(storage.ham_adi[0, :2], direct["H_adi"])
    assert np.allclose(storage.d1ham_adi[0, :2], direct["dH_adi"])
    assert np.allclose(storage.dc1_adi[0, :2], direct["DC1_adi"])
    assert np.allclose(storage.nac_adi[0, :2, 0, 1].real, np.array([1.0, 1.5]) * direct["DC1_adi"][:, 0, 0, 1].real)


def test_shin_metiu_on_the_fly_dvr_returns_adiabatic_model_data():
    model = ShinMetiuDVRModel(params={"N": 41, "r_min": -20.2, "r_max": 19.8, "model": 1})
    result = model.evaluate(np.asarray([[1.3, 2.1]]))

    assert result["H_adi"].shape == (2, 2, 2)
    assert result["dH_adi"].shape == (2, 1, 2, 2)
    assert result["DC1_adi"].shape == (2, 1, 2, 2)
    assert np.allclose(result["H_adi"], np.swapaxes(result["H_adi"].conj(), -1, -2))
    assert np.allclose(result["DC1_adi"][:, 0] + np.swapaxes(result["DC1_adi"][:, 0], -1, -2), 0.0, atol=1.0e-10)


def test_bundled_shin_metiu_dvr_file_loads():
    data = load_shin_metiu_dvr_data(bundled_shin_metiu_dvr_path(1))
    result = data.interpolate(data.R_grid[[10, 20]])

    assert data.nstates == 2
    assert result["H_adi"].shape == (2, 2, 2)
    assert result["mu"].shape == (2, 2, 2)
    assert np.all(np.diff(data.R_grid) > 0.0)


def _synthetic_dvr_data():
    R_grid = np.array([-1.0, 0.0, 1.0])
    eigvals = np.array([[0.0, 0.5], [0.1, 0.6], [0.2, 0.7]])
    d_V = np.zeros((3, 2, 2))
    d_V[:, 0, 0] = 0.1
    d_V[:, 1, 1] = 0.2
    d_V[:, 0, 1] = [0.01, 0.02, 0.03]
    d_V[:, 1, 0] = d_V[:, 0, 1]
    nac = np.zeros((3, 2, 2))
    nac[:, 0, 1] = [0.0, 0.01, 0.02]
    nac[:, 1, 0] = -nac[:, 0, 1]
    mu = np.zeros((3, 2, 2))
    mu[:, 0, 0] = [1.0, 1.2, 1.4]
    mu[:, 1, 1] = [2.0, 2.2, 2.4]
    mu[:, 0, 1] = 0.3
    mu[:, 1, 0] = 0.3
    d_mu = np.zeros_like(mu)
    mu_deri = np.zeros_like(mu)
    mu_deri[:, 0, 0] = 0.2
    mu_deri[:, 1, 1] = 0.4
    mu_deri[:, 0, 1] = 0.05
    mu_deri[:, 1, 0] = 0.05
    return ShinMetiuDVRData(R_grid, eigvals, d_V, nac, mu, d_mu, mu_deri)
