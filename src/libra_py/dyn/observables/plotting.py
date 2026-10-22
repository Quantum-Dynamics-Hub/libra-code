# Copyright (C) 2022-2026 Matthew Dutra, Alexey V. Akimov
# SPDX-License-Identifier: GPL-3.0-or-later

"""Headless Gaussian dynamics plots, with one image file per observable or snapshot.

Plotting is an explicit CPU boundary: NumPy arrays or tensors are accepted,
and tensors are detached only for Matplotlib. Figure objects are closed after
saving. Functions return paths; they never invoke notebook display or show().
"""
from pathlib import Path
import numpy as np
from ..backends.arrays import to_numpy


def _figure():
    # FigureCanvasAgg avoids changing a user's global matplotlib backend.
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig = Figure(figsize=(7, 4.5), layout="constrained")
    FigureCanvasAgg(fig)
    return fig, fig.subplots()


def _write(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    fig.clear()
    return path


def energy_and_pops(result, directory):
    """Save energy.png, energy_change.png, norm.png, populations.png."""
    h = result
    t = to_numpy(h["time"])
    paths = []
    for key, filename, ylabel in (("Etot", "energy", "Energy (Ha)"),
                                  ("dEtot", "energy_change", "E(t) - E(0) (Ha)"),
                                  ("norm", "norm", "Wavefunction norm")):
        fig, ax = _figure()
        ax.plot(t, to_numpy(h[key]))
        ax.set(xlabel="Time (a.u.)", ylabel=ylabel)
        paths.append(_write(fig, Path(directory)/(filename+".png")))
    fig, ax = _figure()
    for k, state in enumerate(to_numpy(h["target_states"][0])):
        ax.plot(t, to_numpy(h["pops"])[:, k], label=f"State {state}")
    ax.set(xlabel="Time (a.u.)", ylabel="Population")
    ax.legend()
    paths.append(_write(fig, Path(directory)/"populations.png"))
    return paths


def trajectories(result, directory, which_states=None):
    """Save separate coordinate, momentum, width and phase plots per DOF."""
    h = result
    labels = to_numpy(h["states"][0])
    states = sorted(set(labels)) if which_states is None else which_states
    paths = []
    for key, name, unit in (("q", "position", "bohr"), ("p", "momentum", "a.u."),
                            ("gaussian_alpha", "width", "bohr^-2"), ("gaussian_phase", "phase", "radian")):
        data = to_numpy(h[key])
        for d in range(data.shape[2]):
            fig, ax = _figure()
            for k, state in enumerate(states):
                idx = np.flatnonzero(labels == state)
                if len(idx):
                    lines = ax.plot(to_numpy(h["time"]), data[:, idx, d], color=f"C{k%10}", alpha=.6)
                    lines[0].set_label(f"State {state}")
            ax.set(xlabel="Time (a.u.)", ylabel=f"{name.capitalize()} DOF {d} ({unit})")
            ax.legend()
            paths.append(_write(fig, Path(directory)/f"{name}_dof{d}.png"))
    return paths


def plot_wf_1D(data, directory):
    """Plot reconstructed density per snapshot, all electronic states together."""
    paths = []
    for k, snap in enumerate(data["snaps"]):
        fig, ax = _figure()
        for j, state in enumerate(data["states"]):
            ax.plot(to_numpy(data["axes"][0]), to_numpy(data["density"][k, j]), label=f"State {state}")
        ax.set(xlabel="Position (bohr)", ylabel="Probability density",
               title=f"t = {float(data['times'][k]):.6g} a.u.")
        ax.legend()
        paths.append(_write(fig, Path(directory)/f"density_step{snap:05d}.png"))
    return paths


def plot_wf_2D(data, directory):
    """Save each snapshot/state density as an individual 2D heatmap."""
    paths = []
    for k, snap in enumerate(data["snaps"]):
        for j, state in enumerate(data["states"]):
            fig, ax = _figure()
            mesh = ax.pcolormesh(to_numpy(data["axes"][0]), to_numpy(data["axes"][1]),
                                 to_numpy(data["density"][k, j]).T, shading="auto")
            fig.colorbar(mesh, ax=ax, label="Probability density")
            ax.set(xlabel="X (bohr)", ylabel="Y (bohr)",
                   title=f"State {state}, t = {float(data['times'][k]):.6g} a.u.")
            paths.append(_write(fig, Path(directory)/f"density_step{snap:05d}_state{state}.png"))
    return paths


def wf_plot(data, directory):
    """Save 1D/2D density images from observables.gaussian.reconstruct_density's returned data."""
    d = len(data["axes"])
    if d not in (1, 2):
        raise ValueError("density plotting supports one or two spatial dimensions")
    return (plot_wf_1D if d == 1 else plot_wf_2D)(data, directory)
