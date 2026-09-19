"""Run the two-state complex-Hermitian gapped CI model."""

import numpy as np
import torch

import libra_py.dynamics.exact_torch.compute as compute


def H_hat_gap(q, params):
    """Return the diabatic Hamiltonian on a two-dimensional grid.

    Parameters
    ----------
    q : torch.Tensor, shape ``(2, nx, ny)``
        Nuclear-coordinate grid.
    params : dict
        Model parameters ``omegax``, ``omegay``, ``kappa``, and ``m``.

    Returns
    -------
    torch.Tensor, shape ``(nx, ny, 2, 2)``
        Complex-Hermitian diabatic Hamiltonian at every grid point.
    """

    x, y = q[0], q[1]
    omegax = params["omegax"]
    omegay = params["omegay"]
    kappa = params["kappa"]
    gap = params["m"]
    common = 0.5 * (omegax**2 * x**2 + omegay**2 * y**2)

    hamiltonian = torch.empty(
        (*x.shape, 2, 2), dtype=torch.complex64, device=x.device
    )
    hamiltonian[..., 0, 0] = common + kappa * x
    hamiltonian[..., 1, 1] = common - kappa * x
    hamiltonian[..., 0, 1] = kappa * y - 1.0j * gap
    hamiltonian[..., 1, 0] = kappa * y + 1.0j * gap
    return hamiltonian


def main():
    shift = 10 ** (-6) * np.sqrt(5 + np.e)
    params = {
        "prefix": "Gapped Model V2",
        "grid_size": [256, 256],
        "q_min": [-3.0 + shift, -3.0 + shift],
        "q_max": [3.0 + shift, 3.0 + shift],
        "save_every_n_steps": 50,
        "dt": 0.05,
        "nsteps": 2400,
        "mass": [1685, 1685],
        "Nstates": 2,
        "representation": "adiabatic",
        "initial_state_index": 0,
        "potential_fn": H_hat_gap,
        "potential_fn_params": {
            "omegax": 0.8,
            "omegay": 0.8,
            "kappa": 1.0,
            "c": 1.0,
            "m": 0.16,
        },
        "psi0_fn": compute.gaussian_wavepacket,
        "psi0_fn_params": {
            "mass": [1685, 1685],
            "alpha": [4, 4],
            "q0": [-1.5, 0.0],
            "p0": [60, 0.0],
        },
        "method": "split-operator",
    }

    solver = compute.exact_tdse_solver_multistate(params)
    solver.solve()


if __name__ == "__main__":
    main()
