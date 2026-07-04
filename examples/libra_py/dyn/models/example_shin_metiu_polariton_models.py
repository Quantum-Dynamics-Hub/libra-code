"""
Calculate and plot Shin-Metiu DVR and polaritonic model quantities.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/models/example_shin_metiu_polariton_models.py
"""

from pathlib import Path
import os

OUT_DIR = Path("shin_metiu_polariton_outputs")
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
from libra_py.dyn.hamiltonians import HamiltonianEngine
from libra_py.dyn.models import (
    ShinMetiuDVRModel,
    ShinMetiuPolaritonModel,
    bundled_shin_metiu_dvr_path,
    load_shin_metiu_dvr_data,
)


def evaluate_model(model, R):
    result = model.evaluate(np.asarray([R]))
    return np.diagonal(result["H_adi"], axis1=-2, axis2=-1).real


def plot_energies(name, R, energies):
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for state in range(energies.shape[-1]):
        ax.plot(R, energies[:, state], label=f"E{state}")
    ax.set_xlabel("R")
    ax.set_ylabel("Energy")
    ax.legend(loc="best")
    fig.tight_layout()
    path = OUT_DIR / f"{name}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def demonstrate_engine(model):
    storage = TensorStorage(
        backend=backend,
        ntraj=1,
        ndof=1,
        nstates=model.nstates,
        ntbf_initial=3,
        ntbf_capacity=3,
    )
    storage.allocate_hamiltonian_derivatives(der_lvl=1)
    storage.q[0, :, 0] = [-1.0, 0.0, 1.0]
    storage.p[0, :, 0] = [5.0, 5.0, 5.0]
    storage.iM[0, :, 0] = 1.0 / 1836.0

    traj = Trajectory(0)
    traj.tbf_ids = [0, 1, 2]
    HamiltonianEngine(backend).evaluate(traj, storage, model, rep="adiabatic")
    return {
        "R": storage.q[0, traj.tbf_ids, 0].copy(),
        "energies": np.diagonal(storage.ham_adi[0, traj.tbf_ids], axis1=-2, axis2=-1).real,
        "dc01": storage.dc1_adi[0, traj.tbf_ids, 0, 0, 1].real,
        "nac01": storage.nac_adi[0, traj.tbf_ids, 0, 1].real,
    }


def main():
    data = load_shin_metiu_dvr_data(bundled_shin_metiu_dvr_path(1))
    R = np.linspace(data.R_grid[200], data.R_grid[-201], 500)

    electronic = ShinMetiuDVRModel(dvr_data=data)
    polariton_2 = ShinMetiuPolaritonModel(
        dvr_data=data,
        model="2-state",
        params={"g_c": 0.005, "omega_c": 0.1, "epsilon": 1.0},
    )
    polariton_4 = ShinMetiuPolaritonModel(
        dvr_data=data,
        model="4-state",
        params={"g_c": 0.005, "omega_c": 0.1, "epsilon": 1.0},
    )

    outputs = []
    for name, model in (
        ("shin_metiu_electronic", electronic),
        ("shin_metiu_polariton_2_state", polariton_2),
        ("shin_metiu_polariton_4_state", polariton_4),
    ):
        energies = evaluate_model(model, R)
        data_path = OUT_DIR / f"{name}.npz"
        np.savez(data_path, R=R, energies=energies)
        plot_path = plot_energies(name, R, energies)
        outputs.append((data_path, plot_path))

    for data_path, plot_path in outputs:
        print(f"saved {data_path} and {plot_path}")

    sample = demonstrate_engine(polariton_4)
    print("\nHamiltonianEngine sample for the 4-state Shin-Metiu polariton")
    print("R:", sample["R"])
    print("energies:\n", sample["energies"])
    print("DC1_01:", sample["dc01"])
    print("velocity-projected NAC01:", sample["nac01"])


if __name__ == "__main__":
    main()
