"""Example 01: plot saved FSSH and Ehrenfest observables."""

import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
PLOT_DIR = HERE / "plots"
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_DIR / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_snapshots(method):
    """Load ordered NPZ snapshots written by FaultTolerantSaver."""

    files = sorted((OUTPUT_DIR / method).glob("step_*.npz"))
    if not files:
        raise FileNotFoundError(
            f"No {method} snapshots found. Run compute.py first."
        )
    snapshots = []
    for path in files:
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
    return np.asarray([snapshot[name] for snapshot in snapshots])


def plot_energies(data):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for axis, (method, snapshots) in zip(axes, data.items()):
        time = series(snapshots, "time")
        for field, label in (
            ("Ekin_ave", "kinetic"),
            ("Epot_ave", "potential"),
            ("Etot_ave", "total"),
        ):
            axis.plot(time, series(snapshots, field), label=label)
        axis.set(title=method.upper(), xlabel="Time (a.u.)")
        axis.grid(alpha=0.25)
        axis.legend()
    axes[0].set_ylabel("Energy (Ha)")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "energies.png", dpi=180)
    plt.close(fig)


def plot_populations(data):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for axis, (method, snapshots) in zip(axes, data.items()):
        time = series(snapshots, "time")
        populations = series(snapshots, "se_pop_adi")
        for state in range(populations.shape[-1]):
            axis.plot(time, populations[:, state], label=f"state {state}")
        axis.set(title=method.upper(), xlabel="Time (a.u.)", ylim=(-0.02, 1.02))
        axis.grid(alpha=0.25)
        axis.legend()
    axes[0].set_ylabel("Adiabatic population")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "populations.png", dpi=180)
    plt.close(fig)


def plot_nuclear_motion(data):
    fig, axes = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    for method, snapshots in data.items():
        time = series(snapshots, "time")
        coordinates = series(snapshots, "q").reshape(len(time), -1)[:, 0]
        momenta = series(snapshots, "p").reshape(len(time), -1)[:, 0]
        axes[0].plot(time, coordinates, label=method)
        axes[1].plot(time, momenta, label=method)
    axes[0].set_ylabel("Coordinate (a.u.)")
    axes[1].set_ylabel("Momentum (a.u.)")
    axes[1].set_xlabel("Time (a.u.)")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "nuclear_motion.png", dpi=180)
    plt.close(fig)


def plot_active_state(snapshots):
    time = series(snapshots, "time")
    states = series(snapshots, "states").reshape(len(time), -1)
    fig, axis = plt.subplots(figsize=(9, 4.5))
    image = axis.imshow(
        states.T,
        aspect="auto",
        interpolation="nearest",
        origin="lower",
        extent=(time[0], time[-1], -0.5, states.shape[1] - 0.5),
        vmin=-0.5,
        vmax=1.5,
        cmap="viridis",
    )
    axis.set(
        xlabel="Time (a.u.)",
        ylabel="Trajectory index",
        title=f"TSH active states ({states.shape[1]} trajectories)",
        yticks=np.arange(states.shape[1]),
    )
    colorbar = fig.colorbar(image, ax=axis, ticks=np.arange(int(states.max()) + 1))
    colorbar.set_label("Active state")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "active_state.png", dpi=180)
    plt.close(fig)


def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    data = {
        "tsh": load_snapshots("tsh"),
        "ehrenfest": load_snapshots("ehrenfest"),
    }
    plot_energies(data)
    plot_populations(data)
    plot_nuclear_motion(data)
    plot_active_state(data["tsh"])
    print("Wrote PNG figures to", PLOT_DIR)


if __name__ == "__main__":
    main()
