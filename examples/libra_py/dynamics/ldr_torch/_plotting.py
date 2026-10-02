"""Shared manuscript-style plots for the three model-specific entrypoints."""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch


def load_pair(folder):
    paths = [folder / "ldr.pt", folder / "dvr.pt"]
    if not all(p.exists() for p in paths):
        raise FileNotFoundError(f"Run this model's run.py first; need {paths[0]} and {paths[1]}")
    pair = [torch.load(p, map_location="cpu", weights_only=True) for p in paths]
    if pair[0]["settings"] != pair[1]["settings"]:
        raise ValueError("LDR and DVR settings differ; rerun both methods with the same settings.py")
    return pair


def curve(ax, pair, values, color, label):
    ax.plot(pair[0]["time"], values[0], color=color, label=label, lw=1.7)
    stride = max(1, len(pair[1]["time"]) // 50)
    ax.plot(pair[1]["time"][::stride], values[1][::stride], "o", color=color,
            markerfacecolor="white", markersize=4, markeredgewidth=0.8)


def finish(fig, axes, path):
    for ax in np.asarray(axes, dtype=object).flat:
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.savefig(path.with_suffix(".png"), dpi=180)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)
    print(f"Saved {path.with_suffix('.png')}")


def plot_tully(pair, folder):
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, layout="constrained")
    fig.suptitle(pair[0]["settings"]["title"] + " — LDR lines, DVR open circles")
    for state, color in enumerate(["tab:orange", "tab:blue"]):
        curve(axes[0], pair, [d["rho"][:, state, state].real for d in pair],
              color, rf"$\rho_{{{state}{state}}}$")
    curve(axes[0], pair, [d["rho"][:, 1, 0].real for d in pair], "tab:purple", r"Re $\rho_{10}$")
    curve(axes[0], pair, [d["rho"][:, 1, 0].imag for d in pair], "tab:green", r"Im $\rho_{10}$")
    axes[0].set_ylabel("Density-matrix elements")
    for key, color in [("kinetic_energy", "tab:purple"), ("potential_energy", "0.3"),
                       ("total_energy", "tab:green")]:
        curve(axes[1], pair, [d[key] for d in pair], color, key.replace("_", " "))
    axes[1].set_ylabel("Energy (Ha)")
    curve(axes[2], pair, [d["position"][:, 0] for d in pair], "tab:cyan", r"$\langle q\rangle$")
    axes[2].set_ylabel("Position (Bohr)")
    axes[2].set_xlabel("Time (a.u.)")
    finish(fig, axes, folder / "comparison")

    count = len(pair[0]["density_time"])
    fig, axes = plt.subplots(1, count, figsize=(4*count, 3.5), squeeze=False, layout="constrained")
    for i, ax in enumerate(axes[0]):
        for data, style, label in zip(pair, ["-", "--"], ["LDR", "DVR"]):
            ax.plot(data["axes"][0], data["density"][i], style, label=label)
        ax.set_title(f"t = {float(pair[0]['density_time'][i]):g} a.u.")
        ax.set_xlabel("q (Bohr)")
    axes[0, 0].set_ylabel("Nuclear probability density")
    finish(fig, axes, folder / "density_snapshots")


def plot_flv(pairs, folder):
    fig, axes = plt.subplots(3, len(pairs), figsize=(7*len(pairs), 10),
                             squeeze=False, layout="constrained")
    fig.suptitle("FLV — LDR lines, DVR open circles")
    for column, (name, pair) in enumerate(pairs.items()):
        axes[0, column].set_title(pair[0]["settings"]["title"])
        for state, color in enumerate(["tab:orange", "tab:blue"]):
            curve(axes[0, column], pair, [d["rho"][:, state, state].real for d in pair],
                  color, rf"$\rho_{{{state}{state}}}$")
        for key, color in [("kinetic_energy", "tab:purple"), ("potential_energy", "0.3"),
                           ("total_energy", "tab:green")]:
            curve(axes[1, column], pair, [d[key] for d in pair], color, key.replace("_", " "))
        for f, color, label in [(0, "tab:cyan", r"$\langle X\rangle$"),
                                (1, "tab:pink", r"$\langle Y\rangle$")]:
            curve(axes[2, column], pair, [d["position"][:, f] for d in pair], color, label)
        axes[0, column].set_ylabel("Population")
        axes[1, column].set_ylabel("Energy (Ha)")
        axes[2, column].set_ylabel("Position (Bohr)")
        axes[2, column].set_xlabel("Time (a.u.)")
        plot_flv_density(pair, folder / name)
    finish(fig, axes, folder / "comparison")


def plot_flv_density(pair, folder):
    count = len(pair[0]["density_time"])
    fig, axes = plt.subplots(2, count, figsize=(4*count, 6), squeeze=False, layout="constrained")
    for col in range(count):
        vmax = max(float(d["density"][col].max()) for d in pair)
        for row, (data, label) in enumerate(zip(pair, ["LDR", "DVR"])):
            x, y = data["axes"]
            ax = axes[row, col]
            mesh = ax.pcolormesh(x, y, data["density"][col].T, shading="auto",
                                 vmin=0, vmax=vmax, cmap="viridis")
            ax.set_title(f"{label}, t={float(data['density_time'][col]):g}")
            ax.set_xlim(*data["settings"]["ldr_bounds"][0])
            ax.set_ylim(*data["settings"]["ldr_bounds"][1])
            ax.set_xlabel("X (Bohr)")
            ax.set_ylabel("Y (Bohr)")
        fig.colorbar(mesh, ax=axes[:, col].tolist(), label="Probability density", shrink=0.85)
    fig.suptitle(pair[0]["settings"]["title"])
    fig.savefig(folder / "density_snapshots.png", dpi=180)
    fig.savefig(folder / "density_snapshots.pdf")
    plt.close(fig)
    print(f"Saved {folder / 'density_snapshots.png'}")


def plot_cases(cases, output_dir):
    parser = argparse.ArgumentParser(description="Plot local LDR/DVR results together")
    parser.add_argument("--quick", action="store_true", help="Read output/quick")
    args = parser.parse_args()
    folder = output_dir / "quick" if args.quick else output_dir
    if "" in cases:
        plot_tully(load_pair(folder), folder)
    else:
        plot_flv({name: load_pair(folder / name) for name in cases}, folder)
