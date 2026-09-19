"""Plot observables loaded from the saved complex-Hermitian TDSE run."""

from pathlib import Path

import matplotlib.pyplot as plt
import torch


DATA_FILE = Path("Gapped Model V2.pt")
PLOT_FILE = Path("Gapped Model V2_observables.png")


def main():
    data = torch.load(DATA_FILE, map_location="cpu", weights_only=False)
    time = data["time"].numpy()
    kinetic = data["kinetic_energy"].numpy()
    potential = data["potential_energy"].numpy()
    total = data["total_energy"].numpy()
    populations_dia = data["rho_dia_all"].diagonal(dim1=-2, dim2=-1).real.numpy()
    populations_adi = data["rho_adi_all"].diagonal(dim1=-2, dim2=-1).real.numpy()

    figure, axes = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True)
    axes[0].plot(time, kinetic, label="Kinetic")
    axes[0].plot(time, potential, label="Potential")
    axes[0].plot(time, total, label="Total", linewidth=2)
    axes[0].set_xlabel("Time (a.u.)")
    axes[0].set_ylabel("Energy (a.u.)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)

    for state in range(populations_dia.shape[1]):
        axes[1].plot(
            time, populations_dia[:, state], "--", label=f"Diabatic {state}"
        )
        axes[1].plot(
            time, populations_adi[:, state], label=f"Adiabatic {state}"
        )
    axes[1].set_xlabel("Time (a.u.)")
    axes[1].set_ylabel("Population")
    axes[1].legend(ncol=2)
    axes[1].grid(alpha=0.25)

    figure.savefig(PLOT_FILE, dpi=200)
    print(f"Saved {PLOT_FILE}")


if __name__ == "__main__":
    main()
