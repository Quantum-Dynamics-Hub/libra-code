from __future__ import annotations

from pathlib import Path
import os

RUN_CACHE = Path(".ehrenfest_tully1_saving_cache")
MPL_CACHE = RUN_CACHE / "matplotlib"
XDG_CACHE = RUN_CACHE / "xdg"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
XDG_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(XDG_CACHE))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from libra_py.dyn.control_params import DynControlParams
from libra_py.dyn.engine import DynamicsEngine
from libra_py.dyn.initialization import (
    make_independent_trajectory_ensemble,
    sample_gaussian_initial_conditions,
)
from libra_py.dyn.models import TullyModel1
from libra_py.dyn.observables import ObservableConfig
from libra_py.dyn.propagation import update_density


def make_output_dir(name):
    out_dir = Path(f"{name}_outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def requested_nsteps(default=1500):
    return int(os.environ.get("EHRENFEST_TULLY1_NSTEPS", str(default)))


def make_tully1_ehrenfest(ntraj=48, seed=7):
    rng = np.random.default_rng(seed)
    q0, p0 = sample_gaussian_initial_conditions(
        ntraj=ntraj,
        q_mean=-9.0,
        q_sigma=0.45,
        p_mean=30.0,
        p_sigma=2.5,
        rng=rng,
    )
    storage, traj = make_independent_trajectory_ensemble(
        q=q0,
        p=p0,
        masses=2000.0,
        amplitudes=[1.0 + 0.0j, 0.0 + 0.0j],
    )
    engine = DynamicsEngine(
        traj,
        storage,
        TullyModel1(),
        params=DynControlParams(force_method=2, ehrenfest_force_option=0),
        method="ehrenfest",
        rep="adiabatic",
    )
    engine.initialize()
    update_density(storage, traj)
    return engine, storage, traj


def full_observable_config():
    return ObservableConfig(
        rep="adiabatic",
        potential="ehrenfest",
        populations=True,
        active_counts=False,
        energies=True,
        include_coordinates=True,
        include_momenta=True,
    )


def compact_observable_config():
    return ObservableConfig(
        rep="adiabatic",
        potential="ehrenfest",
        populations=True,
        active_counts=False,
        energies=True,
        include_coordinates=False,
        include_momenta=False,
    )


def run_with_saver(
    saver,
    config,
    nsteps=None,
    dt=1.0,
    save_stride=1,
    metadata=None,
):
    engine, storage, traj = make_tully1_ehrenfest()
    nsteps = requested_nsteps() if nsteps is None else int(nsteps)
    metadata = metadata or {"model": "TullyModel1", "method": "Ehrenfest"}

    saver.save_observables(storage, traj, config, step=0, time=engine.time, metadata=metadata)
    for step in range(1, nsteps + 1):
        engine.step(dt)
        if step % save_stride == 0:
            saver.save_observables(storage, traj, config, step=step, time=engine.time)
    return nsteps


def add_mean_energies(observables):
    if "kinetic_energy" in observables:
        observables["mean_kinetic_energy"] = np.mean(observables["kinetic_energy"], axis=1)
    if "potential_energy" in observables:
        observables["mean_potential_energy"] = np.mean(observables["potential_energy"], axis=1)
    if "total_energy" in observables:
        observables["mean_total_energy"] = np.mean(observables["total_energy"], axis=1)
    return observables


def save_npz_copy(out_dir, name, observables):
    path = out_dir / name
    np.savez(path, **observables)
    return path


def plot_energy_and_populations(out_dir, observables, prefix="observables"):
    time = observables["time"]
    saved = []

    if "mean_total_energy" in observables:
        fig, ax = plt.subplots(figsize=(7.0, 4.0))
        if "mean_kinetic_energy" in observables:
            ax.plot(time, observables["mean_kinetic_energy"], label="kinetic")
        if "mean_potential_energy" in observables:
            ax.plot(time, observables["mean_potential_energy"], label="Ehrenfest potential")
        ax.plot(time, observables["mean_total_energy"], label="total")
        ax.set_xlabel("time, a.u.")
        ax.set_ylabel("trajectory average, Ha")
        ax.legend(loc="best")
        fig.tight_layout()
        path = out_dir / f"{prefix}_energy.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        saved.append(path)

    if "mean_populations" in observables:
        fig, ax = plt.subplots(figsize=(7.0, 4.0))
        for state in range(observables["mean_populations"].shape[-1]):
            ax.plot(time, observables["mean_populations"][:, state], label=f"state {state}")
        ax.set_xlabel("time, a.u.")
        ax.set_ylabel("population")
        ax.set_ylim(-0.03, 1.03)
        ax.legend(loc="best")
        fig.tight_layout()
        path = out_dir / f"{prefix}_populations.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        saved.append(path)

    return saved


def plot_phase_portrait(out_dir, observables, prefix="observables"):
    if "q" not in observables or "p" not in observables:
        return None
    q = observables["q"][:, :, 0]
    p = observables["p"][:, :, 0]
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    for path_idx in range(q.shape[1]):
        ax.plot(q[:, path_idx], p[:, path_idx], color="0.25", alpha=0.22, linewidth=0.8)
    ax.scatter(q[0], p[0], s=12, color="tab:blue", label="initial")
    ax.scatter(q[-1], p[-1], s=12, color="tab:red", label="final")
    ax.set_xlabel("q")
    ax.set_ylabel("p")
    ax.legend(loc="best")
    fig.tight_layout()
    path = out_dir / f"{prefix}_phase_portraits.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def print_summary(label, out_dir, observables, saved_paths):
    print(label)
    print(f"output directory: {out_dir}")
    for path in saved_paths:
        if path is not None:
            print(f"saved {path}")
    if "mean_populations" in observables:
        print("final mean populations:", observables["mean_populations"][-1])
    if "mean_total_energy" in observables:
        print("initial/final mean total energy:", observables["mean_total_energy"][[0, -1]])
