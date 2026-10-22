"""Run and plot the Tully-1 FSSH/GFSH/FSSH2 decoherence comparison.

By default, simulate all nine cases and write individual PNG figures.
Use --compute-only to skip plotting or --plot-only to reuse saved snapshots.
All physical quantities use atomic units.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from libra_py.dyn.engine import DynamicsEngine, dynamics_defaults
from libra_py.dyn.initialization import (
    make_independent_trajectory_ensemble,
    sample_gaussian_initial_conditions,
)
from libra_py.dyn.models import TullyModel1
from libra_py.dyn.savers.disk import FaultTolerantSaver


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output" / "tsh_decoherence_comparison"
PLOT_DIR = HERE / "plots" / "tsh_decoherence_comparison"
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_DIR / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MASS = 2000.0
DT = 0.5

# C++ tsh_method option numbers are preserved by DynamicsEngine.
TSH_METHODS = {"fssh": 0, "gfsh": 1, "fssh2": 7}
DECOHERENCE = {
    "coherent": {"decoherence_algo": -1, "decoherence_times_type": -1},
    # ID-A acts at every nontrivial attempted hop. Rates are not needed.
    "ida": {
        "decoherence_algo": 1,
        "instantaneous_decoherence_variant": 1,
        "collapse_option": 0,
        "decoherence_times_type": -1,
    },
    # SDM is applied before hop proposal; here its rates come from EDC.
    "sdm": {
        "decoherence_algo": 0,
        "decoherence_times_type": 1,
        "decoherence_C_param": 1.0,
        "decoherence_eps_param": 0.1,
        "sdm_norm_tolerance": 1.0e-10,
    },
}
PROPERTIES = [
    "timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave",
    "dEtot_ave", "states", "se_pop_adi", "sh_pop_adi", "q", "p",
]


def make_system(ntraj):
    """Return the same reproducible incoming ensemble for every method."""

    q, p = sample_gaussian_initial_conditions(
        ntraj=ntraj,
        q_mean=-5.0,
        q_sigma=0.15,
        p_mean=20.0,
        p_sigma=0.5,
        rng=np.random.default_rng(19),
    )
    return make_independent_trajectory_ensemble(
        q=q,
        p=p,
        masses=MASS,
        amplitudes=[1.0 + 0.0j, 0.0 + 0.0j],
        active_state=0,
    )


def run_case(method_name, correction_name, ntraj, nsteps, save_stride):
    """Run and save one TSH/decoherence combination."""

    storage, trajectory = make_system(ntraj)
    options = {
        "tsh_method": TSH_METHODS[method_name],
        "rep_tdse": 1,
        "rep_sh": 1,
        "force_method": 1,
        "hop_acceptance_algo": 20,
        "momenta_rescaling_algo": 200,
        "electronic_integrator": 4,
        "num_electronic_substeps": 4,
        "dt": DT,
        "nsteps": nsteps,
        "nprint": save_stride,
        "properties_to_save": PROPERTIES,
    }
    options.update(DECOHERENCE[correction_name])
    params = dynamics_defaults(options)

    case_name = f"{method_name}_{correction_name}"
    saver = FaultTolerantSaver(OUTPUT_DIR / case_name, mode="w")
    with DynamicsEngine(
        trajectory,
        storage,
        TullyModel1(),
        params=params,
        method="tsh",
        rng=np.random.default_rng(71),
        saver=saver,
        metadata={"tsh_method": method_name, "decoherence": correction_name},
    ) as engine:
        results = engine.run()

    final_states = np.asarray(results[-1].active_states)
    populations = np.bincount(final_states, minlength=2) / ntraj
    print(f"{case_name:14s} final SH populations: {populations}")



COLORS = {"coherent": "black", "ida": "tab:orange", "sdm": "tab:blue"}
LABELS = {"coherent": "coherent", "ida": "ID-A", "sdm": "SDM (EDC)"}


def load_case(method, correction):
    """Load one ordered sequence of fault-tolerant NPZ snapshots."""

    case = f"{method}_{correction}"
    paths = sorted((OUTPUT_DIR / case).glob("step_*.npz"))
    if not paths:
        raise FileNotFoundError(
            f"No snapshots found for {case}. Run "
            "run.py first."
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
    for axis, method in zip(axes, TSH_METHODS):
        for correction in DECOHERENCE:
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
    for axis, method in zip(axes, TSH_METHODS):
        for correction in DECOHERENCE:
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
    for row, method in enumerate(TSH_METHODS):
        for column, correction in enumerate(DECOHERENCE):
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
            if row == len(TSH_METHODS) - 1:
                axis.set_xlabel("Time (a.u.)")
    colorbar = fig.colorbar(image, ax=axes, ticks=[0, 1], shrink=0.85)
    colorbar.set_label("Active state")
    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.07, top=0.95,
                        wspace=0.12, hspace=0.15)
    path = PLOT_DIR / "active_states.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def make_plots():
    """Load all cases and create the comparison PNG files."""

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    data = {
        (method, correction): load_case(method, correction)
        for method in TSH_METHODS
        for correction in DECOHERENCE
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



def parse_args():
    """Read simulation sizes and optionally select only computation or plotting."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ntraj", type=int, default=12)
    parser.add_argument("--nsteps", "--steps", dest="nsteps", type=int, default=1200)
    parser.add_argument("--save-stride", type=int, default=10)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--compute-only", action="store_true",
                      help="save snapshots without making figures")
    mode.add_argument("--plot-only", action="store_true",
                      help="plot existing snapshots without rerunning dynamics")
    args = parser.parse_args()
    if args.ntraj < 1 or args.nsteps < 1 or args.save_stride < 1:
        parser.error("ntraj, nsteps, and save-stride must be positive")
    return args


def main():
    """Simulate all method/correction pairs, then plot their saved observables."""
    args = parse_args()
    if not args.plot_only:
        for method_name in TSH_METHODS:
            for correction_name in DECOHERENCE:
                run_case(method_name, correction_name, args.ntraj,
                         args.nsteps, args.save_stride)
        print("Snapshots written to", OUTPUT_DIR)
    if not args.compute_only:
        make_plots()


if __name__ == "__main__":
    main()
