import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import HamiltonianEngine
from libra_py.dyn.models import (
    ShinMetiuDVRData,
    ShinMetiuPolaritonModel,
    build_four_state_polariton_hamiltonian,
    build_two_state_polariton_hamiltonian,
    dipole_self_energy,
    polariton_info,
)


def test_dipole_self_energy_uses_mu_matrix_product_and_derivative():
    mu = np.array([[1.0, 0.2], [0.2, 2.0]])
    mu_deri = np.array([[0.1, 0.3], [0.3, 0.4]])
    result = dipole_self_energy(mu, mu_deri, g_c=0.02, omega_c=0.1, epsilon=0.5)
    prefactor = 0.5**2 * 0.02**2 / 0.1

    assert np.allclose(result["D_square"], prefactor * (mu @ mu))
    assert np.allclose(result["D_square_deri"], prefactor * (mu @ mu_deri + mu_deri @ mu))


def test_two_state_polariton_builder_matches_documented_block():
    eigvals = np.array([0.1, 0.7])
    H_deri = np.array([[0.2, 0.0], [0.0, 0.3]])
    mu = np.array([[1.0, 0.4], [0.4, 2.0]])
    mu_deri = np.array([[0.1, 0.05], [0.05, 0.2]])
    dse = dipole_self_energy(mu, mu_deri, g_c=0.01, omega_c=0.2, epsilon=1.5)

    H, dH, dc = build_two_state_polariton_hamiltonian(
        eigvals,
        H_deri,
        mu,
        mu_deri,
        dse["D_square"],
        dse["D_square_deri"],
        g_c=0.01,
        omega_c=0.2,
        epsilon=1.5,
    )

    assert np.isclose(H[0, 0], eigvals[1] + dse["D_square"][1, 1])
    assert np.isclose(H[1, 1], eigvals[0] + dse["D_square"][0, 0] + 0.2)
    assert np.isclose(H[0, 1], 0.01 * 1.5 * mu[1, 0])
    assert np.isclose(dH[0, 1], 0.01 * 1.5 * mu_deri[1, 0])
    assert np.allclose(dc, 0.0)


def test_four_state_polariton_builder_is_symmetric_and_has_expected_nac_blocks():
    eigvals = np.array([0.1, 0.7])
    H_deri = np.array([[0.2, 0.01], [0.01, 0.3]])
    nac = np.array([[0.0, 0.04], [-0.04, 0.0]])
    mu = np.array([[1.0, 0.4], [0.4, 2.0]])
    mu_deri = np.array([[0.1, 0.05], [0.05, 0.2]])
    dse = dipole_self_energy(mu, mu_deri, g_c=0.01, omega_c=0.2, epsilon=1.5)

    H, dH, dc = build_four_state_polariton_hamiltonian(
        eigvals,
        H_deri,
        nac,
        mu,
        mu_deri,
        dse["D_square"],
        dse["D_square_deri"],
        g_c=0.01,
        omega_c=0.2,
        epsilon=1.5,
    )

    assert H.shape == (4, 4)
    assert np.allclose(H, H.T)
    assert np.allclose(dc + dc.T, 0.0)
    assert np.isclose(dc[0, 1], nac[0, 1])
    assert np.isclose(dc[2, 3], nac[0, 1])
    assert np.isclose(H[0, 2], 0.01 * 1.5 * mu[0, 0])
    assert np.isclose(H[3, 3], eigvals[1] + dse["D_square"][1, 1] + 1.5 * 0.2)
    assert np.allclose(dH, dH.T)


def test_polariton_info_and_model_return_adiabatic_engine_fields():
    data = _synthetic_dvr_data()
    info = polariton_info(0.0, data, model="4-state", g_c=0.01, omega_c=0.2, epsilon=1.0)

    assert set(info) == {"H_adi", "dH_adi", "DC1_adi"}
    assert info["H_adi"].shape == (4, 4)
    assert info["dH_adi"].shape == (4, 4)
    assert np.allclose(info["H_adi"], info["H_adi"].T)

    model = ShinMetiuPolaritonModel(
        dvr_data=data,
        model="4-state",
        params={"g_c": 0.01, "omega_c": 0.2, "epsilon": 1.0},
    )
    result = model.evaluate(np.asarray([[-0.5, 0.5]]))
    assert result["H_adi"].shape == (2, 4, 4)
    assert result["dH_adi"].shape == (2, 1, 4, 4)
    assert result["DC1_adi"].shape == (2, 1, 4, 4)


def test_polariton_model_writes_velocity_projected_nac_through_engine():
    model = ShinMetiuPolaritonModel(
        dvr_data=_synthetic_dvr_data(),
        model="4-state",
        params={"g_c": 0.01, "omega_c": 0.2},
    )
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=4,
        ntbf_initial=2,
        ntbf_capacity=2,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    storage.q[0, :, 0] = [-0.5, 0.5]
    storage.p[0, :, 0] = [2.0, 4.0]
    storage.iM[0, :, 0] = 0.25
    traj = Trajectory(0)
    traj.tbf_ids = [0, 1]

    HamiltonianEngine(backend).evaluate(traj, storage, model, rep="adiabatic")

    assert np.allclose(storage.ham_adi[0, :2], np.swapaxes(storage.ham_adi[0, :2].conj(), -1, -2))
    expected_velocity = np.array([0.5, 1.0])
    assert np.allclose(storage.nac_adi[0, :2, 0, 1].real, expected_velocity * storage.dc1_adi[0, :2, 0, 0, 1].real)
    assert np.allclose(storage.hvib_adi[0, :2], storage.ham_adi[0, :2] - 1j * storage.nac_adi[0, :2])


def _synthetic_dvr_data():
    R_grid = np.array([-1.0, 0.0, 1.0])
    eigvals = np.array([[0.0, 0.5], [0.1, 0.6], [0.2, 0.7]])
    d_V = np.zeros((3, 2, 2))
    d_V[:, 0, 0] = 0.1
    d_V[:, 1, 1] = 0.2
    d_V[:, 0, 1] = 0.02
    d_V[:, 1, 0] = 0.02
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
