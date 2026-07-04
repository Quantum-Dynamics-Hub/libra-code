"""
Ehrenfest dynamics for a small Tully model 1 trajectory ensemble.

The script is written as consecutive code blocks so the same pieces can be
copied into a notebook.

Run from the repository root:

    PYTHONPATH=src python examples/libra_py/dyn/propagation/example_ehrenfest_tully1_ensemble.py
"""

# %% 1. Imports and output directory
from pathlib import Path
import os

OUT_DIR = Path("ehrenfest_tully1_outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUT_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(OUT_DIR / ".cache"))

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
from libra_py.dyn.savers import HDF5Saver, load_hdf5_steps, stack_saved_observables


# %% 2. Sample initial conditions
rng = np.random.default_rng(7)
ntraj = 48
mass = 2000.0
q0, p0 = sample_gaussian_initial_conditions(
    ntraj=ntraj,
    q_mean=-9.0,
    q_sigma=0.45,
    p_mean=30.0,
    p_sigma=2.5,
    rng=rng,
)


# %% 3. Build storage, model, and Ehrenfest engine
model = TullyModel1()
storage, traj = make_independent_trajectory_ensemble(
    q=q0,
    p=p0,
    masses=mass,
    amplitudes=[1.0 + 0.0j, 0.0 + 0.0j],
)

engine = DynamicsEngine(
    traj,
    storage,
    model,
    params=DynControlParams(force_method=2, ehrenfest_force_option=0),
    method="ehrenfest",
    rep="adiabatic",
)
engine.initialize()
update_density(storage, traj)


# %% 4. Run dynamics and save observables
dt = 1.0
nsteps = int(os.environ.get("EHRENFEST_TULLY1_NSTEPS", "1500"))
observable_config = ObservableConfig(
    rep="adiabatic",
    potential="ehrenfest",
    populations=True,
    active_counts=False,
    energies=True,
    include_coordinates=True,
    include_momenta=True,
)
saver_filename = "observables.hdf"

with HDF5Saver(output_dir=OUT_DIR, filename=saver_filename, mode="w") as saver:
    saver.save_observables(
        storage,
        traj,
        observable_config,
        step=0,
        time=engine.time,
        metadata={"model": "TullyModel1", "method": "Ehrenfest"},
    )
    for step in range(1, nsteps + 1):
        engine.step(dt)
        saver.save_observables(storage, traj, observable_config, step=step, time=engine.time)


# %% 5. Load saved observables for analysis and plotting
records = load_hdf5_steps(output_dir=OUT_DIR, filename=saver_filename)
observables = stack_saved_observables(
    records,
    keys=(
        "time",
        "q",
        "p",
        "populations",
        "mean_populations",
        "kinetic_energy",
        "potential_energy",
        "total_energy",
    ),
)
observables["mean_kinetic_energy"] = np.mean(observables["kinetic_energy"], axis=1)
observables["mean_potential_energy"] = np.mean(observables["potential_energy"], axis=1)
observables["mean_total_energy"] = np.mean(observables["total_energy"], axis=1)
data_path = OUT_DIR / "observables.npz"
np.savez(data_path, **observables)


# %% 6. Plot trajectory-averaged energy
time = observables["time"]
fig, ax = plt.subplots(figsize=(7.0, 4.0))
ax.plot(time, observables["mean_kinetic_energy"], label="kinetic")
ax.plot(time, observables["mean_potential_energy"], label="Ehrenfest potential")
ax.plot(time, observables["mean_total_energy"], label="total")
ax.set_xlabel("time, a.u.")
ax.set_ylabel("trajectory average, Ha")
ax.legend(loc="best")
fig.tight_layout()
energy_plot = OUT_DIR / "energy.png"
fig.savefig(energy_plot, dpi=180)
plt.close(fig)


# %% 7. Plot trajectory-averaged electronic populations
fig, ax = plt.subplots(figsize=(7.0, 4.0))
for state in range(observables["mean_populations"].shape[-1]):
    ax.plot(time, observables["mean_populations"][:, state], label=f"state {state}")
ax.set_xlabel("time, a.u.")
ax.set_ylabel("population")
ax.set_ylim(-0.03, 1.03)
ax.legend(loc="best")
fig.tight_layout()
population_plot = OUT_DIR / "populations.png"
fig.savefig(population_plot, dpi=180)
plt.close(fig)


# %% 8. Plot trajectory-resolved q-p phase portraits
q = observables["q"][:, :, 0]
p = observables["p"][:, :, 0]
fig, ax = plt.subplots(figsize=(6.0, 5.0))
for path in range(ntraj):
    ax.plot(q[:, path], p[:, path], color="0.25", alpha=0.22, linewidth=0.8)
ax.scatter(q[0], p[0], s=12, color="tab:blue", label="initial")
ax.scatter(q[-1], p[-1], s=12, color="tab:red", label="final")
ax.set_xlabel("q")
ax.set_ylabel("p")
ax.legend(loc="best")
fig.tight_layout()
phase_plot = OUT_DIR / "phase_portraits.png"
fig.savefig(phase_plot, dpi=180)
plt.close(fig)


# %% 9. Print a small run summary
print(f"saved {OUT_DIR / saver_filename}")
print(f"saved {data_path}")
print(f"saved {energy_plot}")
print(f"saved {population_plot}")
print(f"saved {phase_plot}")
print("final mean populations:", observables["mean_populations"][-1])
print("initial/final mean total energy:", observables["mean_total_energy"][[0, -1]])
