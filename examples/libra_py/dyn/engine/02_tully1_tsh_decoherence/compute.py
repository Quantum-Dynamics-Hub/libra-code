"""Example 02: compare coherent, ID-A, and SDM TSH on Tully-1."""

import argparse
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
        rng=np.random.default_rng(71),
        saver=saver,
        metadata={"tsh_method": method_name, "decoherence": correction_name},
    ) as engine:
        results = engine.run()

    final_states = np.asarray(results[-1].active_states)
    populations = np.bincount(final_states, minlength=2) / ntraj
    print(f"{case_name:14s} final SH populations: {populations}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ntraj", type=int, default=12)
    parser.add_argument("--nsteps", type=int, default=1200)
    parser.add_argument("--save-stride", type=int, default=10)
    return parser.parse_args()


def main():
    """Run all three hopping methods with coherent, ID-A, and SDM dynamics."""

    args = parse_args()
    for method_name in TSH_METHODS:
        for correction_name in DECOHERENCE:
            run_case(
                method_name,
                correction_name,
                args.ntraj,
                args.nsteps,
                args.save_stride,
            )
    print("Snapshots written to", OUTPUT_DIR)


if __name__ == "__main__":
    main()
