"""Run Tully-1 TSH and Ehrenfest ensembles and save NPZ snapshots."""

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
OUTPUT_DIR = HERE / "output"
NTRAJ = 12
MASS = 2000.0
DT = 0.5
NSTEPS = 1200
SAVE_STRIDE = 5
PROPERTIES = [
    "timestep",
    "time",
    "Ekin_ave",
    "Epot_ave",
    "Etot_ave",
    "dEtot_ave",
    "states",
    "se_pop_adi",
    "sh_pop_adi",
    "q",
    "p",
    "f",
]


def initial_conditions():
    """Use the same reproducible incoming wavepacket for both methods."""

    return sample_gaussian_initial_conditions(
        ntraj=NTRAJ,
        q_mean=-5.0,
        q_sigma=0.15,
        p_mean=20.0,
        p_sigma=0.5,
        rng=np.random.default_rng(19),
    )


def make_system():
    q, p = initial_conditions()
    return make_independent_trajectory_ensemble(
        q=q,
        p=p,
        masses=MASS,
        amplitudes=[1.0 + 0.0j, 0.0 + 0.0j],
        active_state=0,
    )


def run_tsh():
    storage, trajectory = make_system()
    params = dynamics_defaults(
        {
            "tsh_method": 0,             # FSSH
            "rep_tdse": 1,               # adiabatic TDSE amplitudes
            "rep_sh": 1,                 # adiabatic hopping states
            "force_method": 1,           # state-specific forces
            "hop_acceptance_algo": 20,   # require DC rescaling feasibility
            "momenta_rescaling_algo": 200,  # conserve energy along the DC
            "electronic_integrator": 4,  # symmetric two-point Hvib
            "num_electronic_substeps": 4,
            "dt": DT,
            "nsteps": NSTEPS,
            "nprint": SAVE_STRIDE,
            "properties_to_save": PROPERTIES,
        }
    )

    saver = FaultTolerantSaver(output_dir=OUTPUT_DIR / "tsh", mode="w")
    with DynamicsEngine(
        trajectory,
        storage,
        TullyModel1(),
        params=params,
        rng=np.random.default_rng(7),
        saver=saver,
    ) as engine:
        results = engine.run()
    print("TSH final states:", results[-1].accepted_states.tolist())
    print("TSH snapshots:", saver.output_dir)


def run_ehrenfest():
    storage, trajectory = make_system()
    params = dynamics_defaults(
        {
            "tsh_method": -1,
            "force_method": 2,
            "rep_tdse": 1,
            "electronic_integrator": 4,
            "num_electronic_substeps": 4,
            "dt": DT,
            "nsteps": NSTEPS,
            "nprint": SAVE_STRIDE,
            "properties_to_save": PROPERTIES,
        }
    )
    saver = FaultTolerantSaver(output_dir=OUTPUT_DIR / "ehrenfest", mode="w")
    with DynamicsEngine(
        trajectory,
        storage,
        TullyModel1(),
        params=params,
        method="ehrenfest",
        saver=saver,
    ) as engine:
        results = engine.run()
    print("Ehrenfest final force:", results[-1].forces)
    print("Ehrenfest snapshots:", saver.output_dir)


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    run_tsh()
    run_ehrenfest()
    print("Run `python plot_dynamics.py` to create PNG figures.")
