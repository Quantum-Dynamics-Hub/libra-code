"""
Calculate and plot Tully model Hamiltonian quantities.

This example exercises the new analytical models in two ways:

* exact-dynamics ordering: Q has shape (ndof, *grid) and the diabatic
  potential has shape (*grid, nstates, nstates)
* trajectory dynamics ordering: the same model object is passed through the
  Hamiltonian machinery to derive adiabatic energies and derivative couplings

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/models/example_tully_models.py
"""

from pathlib import Path
import os

OUT_DIR = Path("tully_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from libra_py.dyn.backends import backend
from libra_py.dyn.core.storage import TensorStorage
from libra_py.dyn.core.trajectory import Trajectory
from libra_py.dyn.hamiltonians import (
    HamiltonianEngine,
    compute_adiabatic_from_diabatic,
    compute_diabatic,
)
from libra_py.dyn.models import TullyModel1, TullyModel2, TullyModel3


def compute_curves(model, x):
    result = compute_hamiltonian_curves(model, x)
    return {
        "x": x,
        "dia_0": result["H_dia"][:, 0, 0].real,
        "dia_1": result["H_dia"][:, 1, 1].real,
        "dia_coupling": result["H_dia"][:, 0, 1].real,
        "adi_0": result["H_adi"][:, 0, 0].real,
        "adi_1": result["H_adi"][:, 1, 1].real,
        "dc_01": result["DC1_adi"][:, 0, 0, 1].real,
    }


def compute_hamiltonian_curves(model, x):
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=len(x),
        ntbf_capacity=len(x),
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)

    traj = Trajectory(0)
    traj.tbf_ids = list(range(len(x)))
    storage.q[0, traj.tbf_ids, 0] = x

    compute_diabatic(storage, traj, model, der_lvl=1)
    compute_adiabatic_from_diabatic(storage, traj, der_lvl=1)
    return {
        "H_dia": storage.ham_dia[0, traj.tbf_ids].copy(),
        "H_adi": storage.ham_adi[0, traj.tbf_ids].copy(),
        "DC1_adi": storage.dc1_adi[0, traj.tbf_ids].copy(),
    }


def plot_curves(name, curves):
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.0), sharex=True)

    axes[0].plot(curves["x"], curves["dia_0"], label="Hdia 00")
    axes[0].plot(curves["x"], curves["dia_1"], label="Hdia 11")
    axes[0].plot(curves["x"], curves["dia_coupling"], label="Hdia 01")
    axes[0].set_ylabel("Diabatic energy")
    axes[0].legend(loc="best")

    axes[1].plot(curves["x"], curves["adi_0"], label="E0")
    axes[1].plot(curves["x"], curves["adi_1"], label="E1")
    axes[1].set_ylabel("Adiabatic energy")
    axes[1].legend(loc="best")

    axes[2].plot(curves["x"], curves["dc_01"], label="d01")
    axes[2].axhline(0.0, color="0.65", linewidth=0.8)
    axes[2].set_xlabel("Nuclear coordinate")
    axes[2].set_ylabel("Derivative coupling")
    axes[2].legend(loc="best")

    fig.suptitle(name)
    fig.tight_layout()
    path = OUT_DIR / f"{name}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def save_curves(name, curves):
    path = OUT_DIR / f"{name}.npz"
    np.savez(path, **curves)
    return path


def demonstrate_hamiltonian_engine(model):
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=2,
        ntbf_initial=3,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)

    traj = Trajectory(0)
    traj.tbf_ids = [0, 1, 2]
    storage.q[0, traj.tbf_ids, 0] = [-6.0, 0.0, 6.0]
    storage.p[0, traj.tbf_ids, 0] = [12.0, 12.0, 12.0]
    storage.iM[0, traj.tbf_ids, 0] = 1.0 / 2000.0

    engine = HamiltonianEngine(backend)
    engine.evaluate(traj, storage, model, rep="adiabatic")
    return {
        "q": storage.q[0, traj.tbf_ids, 0].copy(),
        "H_adi": storage.ham_adi[0, traj.tbf_ids].copy(),
        "DC1_adi": storage.dc1_adi[0, traj.tbf_ids].copy(),
        "NAC_adi": storage.nac_adi[0, traj.tbf_ids].copy(),
        "Hvib_adi": storage.hvib_adi[0, traj.tbf_ids].copy(),
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    x = np.linspace(-12.0, 12.0, 1201)
    models = {
        "tully_model1": TullyModel1(),
        "tully_model2": TullyModel2(),
        "tully_model3": TullyModel3(),
    }

    for name, model in models.items():
        curves = compute_curves(model, x)
        data_path = save_curves(name, curves)
        plot_path = plot_curves(name, curves)
        print(f"{name}: saved {data_path} and {plot_path}")

    sample = demonstrate_hamiltonian_engine(models["tully_model1"])
    print("\nHamiltonianEngine sample for Tully model 1")
    print("q:", sample["q"])
    print("adiabatic energies:\n", np.diagonal(sample["H_adi"], axis1=-2, axis2=-1).real)
    print("adiabatic d01:\n", sample["DC1_adi"][:, 0, 0, 1].real)
    print("time-derivative NAC01:\n", sample["NAC_adi"][:, 0, 1].real)


if __name__ == "__main__":
    main()
