"""Run a reproducible Tully-1 FSSH ensemble through DynamicsEngine.

All quantities use atomic units. Run with --steps 10 for a quick check.
The saved populations are electronic (adiabatic) ensemble averages.
"""
import argparse
from pathlib import Path
import numpy as np

from libra_py.dyn.engine import DynamicsEngine, dynamics_defaults
from libra_py.dyn.initialization import (
    make_independent_trajectory_ensemble, sample_gaussian_initial_conditions,
)
from libra_py.dyn.models import TullyModel1
from libra_py.dyn.savers import (
    FaultTolerantSaver, HDF5TimeSeriesSaver, load_saved_steps,
    load_hdf5_timeseries, stack_saved_observables,
)
from libra_py.dyn.observables import plotting

METHOD = 'tsh'
HERE = Path(__file__).resolve().parent


def configuration(steps=1200, save_stride=5):
    """Specify integration and output options; initial conditions are below."""
    return dynamics_defaults(dict(
        dt=0.5, nsteps=steps, nprint=save_stride,
        tsh_method=0, force_method=1,
        rep_tdse=1, rep_sh=1,  # adiabatic amplitudes and hopping labels
        electronic_integrator=4, num_electronic_substeps=4,
        hop_acceptance_algo=20, momenta_rescaling_algo=200,
        # Hops conserve energy by rescaling along the derivative coupling.
        properties_to_save=["timestep", "time", "Ekin_ave", "Epot_ave",
                            "Etot_ave", "dEtot_ave", "se_pop_adi",
                            "sh_pop_adi", "states", "q", "p", "f"],
    ))


def make_system(ntraj=12):
    """Sample the same incoming packet for FSSH and Ehrenfest comparisons."""
    q, p = sample_gaussian_initial_conditions(
        ntraj=ntraj,
        q_mean=-5.0, q_sigma=0.15,  # incoming position and its standard deviation
        p_mean=20.0, p_sigma=0.5,   # incoming momentum and its standard deviation
        rng=np.random.default_rng(19),
    )
    return make_independent_trajectory_ensemble(
        q=q, p=p, masses=2000.0,
        amplitudes=[1.0+0.0j, 0.0+0.0j],  # initially lower adiabatic state
        active_state=0,
    )


def make_plots(history, output):
    """Use the same filenames as QTAG, retaining ensemble-specific semantics."""
    populations = history["se_pop_adi"]
    summary = dict(time=history["time"], Etot=history["Etot_ave"],
                   dEtot=history["Etot_ave"]-history["Etot_ave"][0],
                   norm=populations.sum(axis=-1), pops=populations,
                   target_states=[[0, 1]])
    # dEtot_ave is an energy spread, not drift; subtract E(0) explicitly.
    plotting.energy_and_pops(summary, output)
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fields = [("q", "position_dof0", "Position (bohr)"),
              ("p", "momentum_dof0", "Momentum (a.u.)")]
    if METHOD == "tsh":
        fields += [("sh_pop_adi", "surface_hopping_populations", "Active-state population"),
                   ("states", "active_states", "Active state")]
    for key, filename, label in fields:
        fig = Figure(figsize=(7, 4.5), layout="constrained")
        FigureCanvasAgg(fig)
        ax = fig.subplots()
        values = np.asarray(history[key]).reshape(len(history["time"]), -1)
        ax.plot(history["time"], values, alpha=0.6)
        ax.set(xlabel="Time (a.u.)", ylabel=label)
        fig.savefig(output / (filename+".png"), dpi=150)
        fig.clear()
    return summary


def main(steps=1200, ntraj=12, save_stride=5, output=None, save_format="npz"):
    """Initialize, propagate, reload saved observables, and write individual PNGs."""
    if steps < 0 or ntraj < 1 or save_stride < 1:
        raise ValueError("steps must be nonnegative; ntraj and save_stride must be positive")
    output = Path(output) if output else HERE / "output" / "tully1_numpy"
    params = configuration(steps, save_stride)
    storage, trajectory = make_system(ntraj)
    if save_format == "hdf5":
        saver = HDF5TimeSeriesSaver(output/"trajectory.hdf", mode="w", flush_stride=100)
    elif save_format == "npz":
        saver = FaultTolerantSaver(output/"snapshots", mode="w", compressed=True)
    else:
        raise ValueError("save_format must be npz or hdf5")
    with DynamicsEngine(trajectory, storage, TullyModel1(), method=METHOD,
                        params=params, rng=np.random.default_rng(7), saver=saver) as engine:
        result = engine.run()
        # Include the final physical time even when it falls between save strides.
        if steps % save_stride:
            saver.save_observables(storage, trajectory, engine.observable_config,
                                   step=storage.timestep, time=engine.time)
    records = (load_saved_steps(output/"snapshots") if save_format == "npz"
               else load_hdf5_timeseries(output/"trajectory.hdf"))
    history = stack_saved_observables(records) if save_format == "npz" else records
    summary = make_plots(history, output)
    print(f"Method: {METHOD}; backend: numpy; model: tully1; completed steps: {steps}")
    print("Final populations:", summary["pops"][-1])
    print("Final norm:", summary["norm"][-1])
    print("Energy change (Ha):", summary["dEtot"][-1])
    print("Saved:", output)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--ntraj", type=int, default=12)
    parser.add_argument("--save-stride", type=int, default=5)
    parser.add_argument("--output")
    parser.add_argument("--save-format", choices=["npz", "hdf5"], default="npz")
    main(**vars(parser.parse_args()))
