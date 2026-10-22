# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""QTAG basics tutorial as a script: NumPy/PyTorch, no notebook or C++ types.

From the repository root:
    PYTHONPATH=src python examples/libra_py/dyn/methods/qtag/1_basics/run.py
    PYTHONPATH=src python examples/libra_py/dyn/methods/qtag/1_basics/run.py --backend torch

Use --steps 10 for a quick smoke run or --case ferretti for the 2D example.
The complete defaults retain the notebook's 800 (1D) or 100 (2D) timesteps.
Every plot is saved to an individual PNG; numerical data is saved in NPZ.
"""
import argparse
from contextlib import nullcontext
from pathlib import Path
import numpy as np
from libra_py.dyn.methods.qtag import initialize, compute, save, plot
from libra_py.dyn.methods.qtag._backend import to_numpy
from models import compute_model


def main(backend="numpy", case="holstein", steps=None, output=None, device="cpu", method=None):
    """Run initialization, dynamics, wavefunction reconstruction, and plotting."""
    output = Path(output or f"qtag_{case}_{backend}")
    # The notebook's two 1D choices (custom/built-in) describe the same model.
    # One pure-array callback covers both here. Its parameters are unchanged.
    if case == "holstein":
        model_params = {
            "case": case, "backend": backend,
            "E_n": [0.0, 15.8114], "x_n": [0.0, 1.0], "k_n": [10.0, 10.0],
            "V": [[0.0, 1.0], [1.0, 0.0]],
            "alpha": [[0.0, 1.5811], [1.5811, 0.0]],
            "x_nm": [[0.0, 2.0], [2.0, 0.0]],
            "ex_cpl_A": 1.0, "ex_cpl_B": 1.5811, "ex_cpl_C": 2.0,
        }
        # Initial packet widths are GWP widths; the BASIS uses QTAG widths.
        # A basis width of 6 means exp(-3*(r-q)^2), whereas wfc_a0=1 means
        # the initial state is proportional to exp(-(r-q0)^2).
        init_params = {
            "wfc_q0": [-2.0], "wfc_p0": [0.0], "wfc_a0": [1.0], "wfc_s0": [0.0],
            "grid_dims": [25], "alp_scl": [6.0], "rho_cut": 1e-12,
        }
        masses = [1.0]
        dt, default_steps, approximation = 0.01, 800, 3  # BAT + exact Gaussian coupling
        grid = dict(xmin=[-4.], xmax=[4.], npoints=[100])
    else:
        model_params = {
            "case": case, "backend": backend,
            "X1": 4.0, "X2": 3.0, "X3": 3.0, "Kx": 0.02, "Ky": 0.01,
            "Delta": 0.01, "alpha": 3.0, "beta": 1.5, "gamma": 0.005,
        }
        # Ky=0.01 matches the notebook's executable cell (its prose says 0.10).
        init_params = {
            "wfc_q0": [0.0, 2.0], "wfc_p0": [0.0, 0.0],
            "wfc_a0": [1.0, 1.0], "wfc_s0": [0.0, 0.0],
            "grid_dims": [5, 5], "alp_scl": [6.0, 6.0], "rho_cut": 1e-12,
        }
        masses = [20000.0, 6667.0]
        dt, default_steps, approximation = 1.0, 100, 0  # BAT
        grid = dict(xmin=[-2., -4.], xmax=[6., 4.], npoints=[25, 25])
    init_params.update(states=[0, 1], init_placement=0, backend=backend, device=device)
    nsteps = default_steps if steps is None else steps
    if method is not None:
        approximation = method
    if case == "ferretti" and approximation in (2, 3):
        raise ValueError("Ferretti coupling is not the separable Gaussian used by exact-coupling methods")

    # No gradients are needed for plotting a long trajectory. Scalar autograd
    # remains available in core/compute when called without this context.
    context = nullcontext()
    if backend == "torch":
        import torch
        inv_mass = torch.tensor([1/m for m in masses], dtype=torch.float64, device=device)
        context = torch.no_grad()
    else:
        inv_mass = np.array([1/m for m in masses], dtype=float)
    with context:
        ntraj, q, p, a, s, labels = initialize.gbfs(init_params)
        q0, p0, a0, s0 = initialize.wfc_centers(init_params)
        coeff = initialize.coeffs(q0, p0, a0, s0, [0], q, p, a, s, labels)

        print(coeff)

        # Same basis updates and synchronization choices as the tutorial.
        # A stable spectral cutoff replaces inversion of nearly singular S.
        dyn_params = {
            "dt": dt, "nsteps": nsteps, "states": [0, 1], "target_states": [0, 1],
            "qtag_pot_approx_method": approximation,
            "q_update_method": 1, "p_update_method": 1, "a_update_method": 1,
            "s_update_method": 0, "q_sync_method": 1, "p_sync_method": 0,
            "a_sync_method": 0, "s_sync_method": 0,
            "mom_calc_type": 1, "d_weight": 1, "decpl_den": 0.01,
            # Match the notebook's iterative momentum fit and its tolerance.
            "momentum_fit": "legacy", "linfit_beta": 1e-3,
            "overlap_rcond": 1e-12,
        }
        result = compute.run_qtag(q, p, a, s, labels, coeff, inv_mass,
                                  dyn_params, compute_model, model_params)
        archive = save.save_qtag_data(result, output/"trajectory.npz")
        save.save_text(result, output)
        plot.energy_and_pops(result, output)
        plot.trajectories(result, output)

        # Tutorial rows [0,199,399,599,799] contain completed steps
        # [1,200,400,600,800]: its files omit the initial state, ours retain it.
        # Match those physical times, also supporting shorter runs and 2D.
        snaps = sorted(set([1] + [int(nsteps*f) for f in (.25, .5, .75, 1)])) if nsteps else [0]
        snaps = [k for k in snaps if k > 0] if nsteps else snaps
        grid.update(snaps=snaps, which_states=[0, 1])
        wave = compute.wf_calc_nD(result, grid)
        plot.wf_plot(wave, output)
        # Plot diabatic and adiabatic potential cuts, one file per spatial DOF.
        for d in range(len(masses)):
            x = np.linspace(grid["xmin"][d], grid["xmax"][d], 200)
            points = np.zeros((len(masses), len(x)))
            points[d] = x
            if backend == "torch":
                points = torch.tensor(points, dtype=torch.float64, device=device)
            h = to_numpy(compute_model(points, model_params).H)
            fig, ax = plot._figure()
            for state in range(2):
                ax.plot(x, h[:, state, state].real, label=f"Diabatic {state}")
                ax.plot(x, np.linalg.eigvalsh(h)[:, state], "--", label=f"Adiabatic {state}")
            ax.set(xlabel=f"DOF {d} (bohr)", ylabel="Potential energy (Ha)")
            ax.legend()
            plot._write(fig, output/f"potential_dof{d}.png")
    print(f"Backend: {backend}; model: {case}; packets: {ntraj}; completed steps: {nsteps}")
    print("Final populations:", to_numpy(result.history["pops"][-1]))
    print("Final norm:", float(result.history["norm"][-1]))
    print("Energy change (Ha):", float(result.history["dEtot"][-1]))
    print("Saved:", archive)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=["numpy", "torch"], default="numpy")
    parser.add_argument("--case", choices=["holstein", "ferretti"], default="holstein")
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--output", default=None)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--method", type=int, choices=range(4), default=None)
    main(**vars(parser.parse_args()))
