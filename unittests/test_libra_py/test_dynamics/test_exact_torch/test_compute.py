from pathlib import Path
import sys

import numpy as np
import pytest
import torch


ROOT = Path(__file__).resolve().parents[4]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from libra_py.dynamics.exact_torch.compute import exact_tdse_solver_multistate


def gapped_complex_hermitian_model(q, params):
    """Complex-Hermitian two-state model used by the regression tests."""

    x, y = q[0], q[1]
    omega_x = params["omegax"]
    omega_y = params["omegay"]
    kappa = params["kappa"]
    gap = params["m"]
    common = 0.5 * (omega_x**2 * x**2 + omega_y**2 * y**2)

    hamiltonian = torch.empty(
        (*x.shape, 2, 2), dtype=torch.complex64, device=x.device
    )
    hamiltonian[..., 0, 0] = common + kappa * x
    hamiltonian[..., 1, 1] = common - kappa * x
    hamiltonian[..., 0, 1] = kappa * y - 1.0j * gap
    hamiltonian[..., 1, 0] = kappa * y + 1.0j * gap
    return hamiltonian


def gaussian(q, params):
    displacement = q - torch.tensor(params["q0"]).view(2, 1, 1)
    phase = 1.0j * torch.sum(
        torch.tensor(params["p0"]).view(2, 1, 1) * q, dim=0
    )
    wavefunction = torch.exp(-torch.sum(displacement**2, dim=0) + phase)
    norm = torch.sqrt(torch.sum(torch.abs(wavefunction) ** 2) * params["dV"])
    return wavefunction / norm


def make_solver():
    grid_size = [12, 10]
    q_min = [-3.0, -3.0]
    q_max = [3.0, 3.0]
    dq = [(hi - lo) / (n - 1) for lo, hi, n in zip(q_min, q_max, grid_size)]
    return exact_tdse_solver_multistate(
        {
            "grid_size": grid_size,
            "q_min": q_min,
            "q_max": q_max,
            "save_every_n_steps": 5,
            "dt": 0.002,
            "nsteps": 20,
            "mass": [1685.0, 1685.0],
            "Nstates": 2,
            "representation": "adiabatic",
            "initial_state_index": 0,
            "potential_fn": gapped_complex_hermitian_model,
            "potential_fn_params": {
                "omegax": 0.8,
                "omegay": 0.8,
                "kappa": 1.0,
                "m": 0.16,
            },
            "psi0_fn": gaussian,
            "psi0_fn_params": {
                "q0": [-1.5, 0.0],
                "p0": [60.0, 0.0],
                "dV": np.prod(dq),
            },
            "method": "split-operator",
            "device": torch.device("cpu"),
        }
    )


def test_complex_hermitian_potential_exponential_and_basis_round_trip():
    solver = make_solver()
    solver.initialize_grids()
    solver.initialize_operators()

    np.testing.assert_allclose(
        solver.expV_half.numpy(),
        torch.linalg.matrix_exp(-0.5j * solver.dt * solver.V).numpy(),
        rtol=2.0e-5,
        atol=2.0e-6,
    )

    adiabatic = solver.psi_r_adi.clone()
    solver.update_dia_r()
    solver.update_adi_r()
    np.testing.assert_allclose(
        solver.psi_r_adi.numpy(), adiabatic.numpy(), rtol=2.0e-5, atol=2.0e-6
    )


def test_complex_hermitian_split_operator_conserves_norm_and_energy():
    solver = make_solver()
    solver.initialize_grids()
    solver.initialize_operators()
    solver.propagate()

    norms = solver.norm.numpy()
    energies = solver.total_energy.numpy()
    np.testing.assert_allclose(norms, norms[0], rtol=2.0e-5, atol=2.0e-6)
    assert np.ptp(energies) < 2.0e-4


def test_nonhermitian_potential_is_rejected():
    solver = make_solver()

    def nonhermitian(q, params):
        value = gapped_complex_hermitian_model(q, params)
        value[..., 0, 1] += 0.1j
        return value

    solver.potential_fn = nonhermitian
    solver.initialize_grids()

    with pytest.raises(ValueError, match="Hermitian"):
        solver.initialize_operators()
