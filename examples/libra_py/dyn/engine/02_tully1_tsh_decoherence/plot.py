"""Example 02: plot the FSSH/GFSH/FSSH2 decoherence comparison."""

import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output" / "tsh_decoherence_comparison"
PLOT_DIR = HERE / "plots" / "tsh_decoherence_comparison"
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_DIR / ".matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


METHODS = ("fssh", "gfsh", "fssh2")
CORRECTIONS = ("coherent", "ida", "sdm")
COLORS = {"coherent": "black", "ida": "tab:orange", "sdm": "tab:blue"}
LABELS = {"coherent": "coherent", "ida": "ID-A", "sdm": "SDM (EDC)"}


def load_case(method, correction):
    """Load one ordered sequence of fault-tolerant NPZ snapshots."""

    case = f"{method}_{correction}"
    paths = sorted((OUTPUT_DIR / case).glob("step_*.npz"))
    if not paths:
        raise FileNotFoundError(
            f"No snapshots found for {case}. Run "
            "compute.py first."
        )
    snapshots = []
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            snapshot = {
                key: np.array(data[key], copy=True)
                for key in data.files
                if key != "__scalars__"
            }
            snapshot.update(json.loads(str(data["__scalars__"])))
            snapshots.append(snapshot)
    return snapshots


def series(snapshots, name):
    """Stack one saved observable along the time axis."""

    return np.asarray([snapshot[name] for snapshot in snapshots])


def state_one_population(values):
    """Extract state-1 population from a saved population time series."""

    values = np.asarray(values)
    if values.shape[-1] < 2:
        raise ValueError("the Tully-1 comparison requires two electronic states")
    return values[..., 1]


def plot_populations(data, observable, filename, ylabel):
    """Plot one population definition in method-resolved panels."""

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8), sharex=True, sharey=True)
    for axis, method in zip(axes, METHODS):
        for correction in CORRECTIONS:
            snapshots = data[(method, correction)]
            axis.plot(
                series(snapshots, "time"),
                state_one_population(series(snapshots, observable)),
                color=COLORS[correction],
                label=LABELS[correction],
            )
        axis.set_title(method.upper())
        axis.set_xlabel("Time (a.u.)")
        axis.set_ylim(-0.02, 1.02)
        axis.grid(alpha=0.25)
    axes[0].set_ylabel(ylabel)
    axes[-1].legend(loc="best")
    fig.tight_layout()
    path = PLOT_DIR / filename
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_energy_drift(data):
    """Compare mean total-energy drift for all nine calculations."""

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8), sharex=True, sharey=True)
    for axis, method in zip(axes, METHODS):
        for correction in CORRECTIONS:
            snapshots = data[(method, correction)]
            time = series(snapshots, "time")
            # dEtot_ave is the ensemble standard deviation, not a difference
            # from the initial energy. Build the actual conservation error.
            energy = series(snapshots, "Etot_ave")
            drift = energy - energy[0]
            axis.plot(time, drift, color=COLORS[correction], label=LABELS[correction])
        axis.set_title(method.upper())
        axis.set_xlabel("Time (a.u.)")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel(r"$\langle E_{tot}(t)-E_{tot}(0)\rangle$ (Ha)")
    axes[-1].legend(loc="best")
    fig.tight_layout()
    path = PLOT_DIR / "energy_drift.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_active_states(data):
    """Plot active-state histories for every trajectory and calculation."""

    fig, axes = plt.subplots(3, 3, figsize=(13.0, 9.0), sharex=True, sharey=True)
    image = None
    for row, method in enumerate(METHODS):
        for column, correction in enumerate(CORRECTIONS):
            axis = axes[row, column]
            snapshots = data[(method, correction)]
            time = series(snapshots, "time")
            states = np.asarray(series(snapshots, "states"))
            states = np.squeeze(states)
            if states.ndim == 1:
                states = states[:, None]
            image = axis.imshow(
                states.T,
                origin="lower",
                aspect="auto",
                interpolation="nearest",
                extent=(time[0], time[-1], -0.5, states.shape[1] - 0.5),
                vmin=0,
                vmax=1,
                cmap="viridis",
            )
            if row == 0:
                axis.set_title(LABELS[correction])
            if column == 0:
                axis.set_ylabel(f"{method.upper()}\ntrajectory")
            if row == len(METHODS) - 1:
                axis.set_xlabel("Time (a.u.)")
    colorbar = fig.colorbar(image, ax=axes, ticks=[0, 1], shrink=0.85)
    colorbar.set_label("Active state")
    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.07, top=0.95,
                        wspace=0.12, hspace=0.15)
    path = PLOT_DIR / "active_states.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main():
    """Load all cases and create the comparison PNG files."""

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    data = {
        (method, correction): load_case(method, correction)
        for method in METHODS
        for correction in CORRECTIONS
    }
    paths = [
        plot_populations(data, "sh_pop_adi", "surface_hopping_population.png",
                         "State-1 SH population"),
        plot_populations(data, "se_pop_adi", "electronic_population.png",
                         "State-1 electronic population"),
        plot_energy_drift(data),
        plot_active_states(data),
    ]
    for path in paths:
        print("Saved", path)


if __name__ == "__main__":
    main()
