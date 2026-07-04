from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import compute_adiabatic_from_diabatic, compute_diabatic


def compute_curves(model, q_values):
    q_values = np.asarray(q_values, dtype=float)
    if q_values.ndim == 1:
        q_values = q_values.reshape(1, -1)
    ndof, npoints = q_values.shape

    native = model.evaluate(q_values)
    nstates = native["H_dia"].shape[-1]
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=ndof,
        nstates=nstates,
        ntbf_initial=npoints,
        ntbf_capacity=npoints,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)

    traj = Trajectory(0)
    traj.tbf_ids = list(range(npoints))
    storage.q[0, traj.tbf_ids, :] = q_values.T

    compute_diabatic(storage, traj, model, der_lvl=1)
    compute_adiabatic_from_diabatic(storage, traj, der_lvl=1)
    return {
        "H_dia": storage.ham_dia[0, traj.tbf_ids].copy(),
        "H_adi": storage.ham_adi[0, traj.tbf_ids].copy(),
        "DC1_adi": storage.dc1_adi[0, traj.tbf_ids].copy(),
    }


def plot_model_curves(name, q_axis, curves, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    q_axis = np.asarray(q_axis)
    H_dia = curves["H_dia"]
    H_adi = curves["H_adi"]
    DC1_adi = curves["DC1_adi"]
    nstates = H_dia.shape[-1]

    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.0), sharex=True)
    for state in range(nstates):
        axes[0].plot(q_axis, H_dia[:, state, state].real, label=f"Hdia {state}{state}")
    if nstates > 1:
        axes[0].plot(q_axis, H_dia[:, 0, 1].real, "--", label="Hdia 01")
    axes[0].set_ylabel("Diabatic energy")
    axes[0].legend(loc="best")

    for state in range(nstates):
        axes[1].plot(q_axis, H_adi[:, state, state].real, label=f"E{state}")
    axes[1].set_ylabel("Adiabatic energy")
    axes[1].legend(loc="best")

    if nstates > 1:
        axes[2].plot(q_axis, DC1_adi[:, 0, 0, 1].real, label="d01")
    else:
        axes[2].plot(q_axis, np.zeros_like(q_axis), label="d00")
    axes[2].axhline(0.0, color="0.65", linewidth=0.8)
    axes[2].set_xlabel("Nuclear coordinate")
    axes[2].set_ylabel("Derivative coupling")
    axes[2].legend(loc="best")

    fig.suptitle(name)
    fig.tight_layout()
    plot_path = out_dir / f"{name}.png"
    data_path = out_dir / f"{name}.npz"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    np.savez(data_path, q=q_axis, **curves)
    return data_path, plot_path
