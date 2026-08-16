"""Examples of Wigner sampling from normal-mode files and a Hessian.

Run from the repository root with::

    PYTHONPATH=src python examples/libra_py/examples_wigner.py
"""

from pathlib import Path

import numpy as np

from libra_py import units
from libra_py.wigner import (
    build_modes_from_hessian,
    generate_wigner_from_hessian,
    prepare_wigner_from_modes,
)


REFERENCE_DIR = Path(__file__).with_name("wigner_reference")
MASSES_DA = {"O": 15.999, "H": 1.00784}


def read_xyz_geometry(filename):
    """Return labels and flattened coordinates in bohr from an XYZ file."""
    lines = Path(filename).read_text().splitlines()
    natoms = int(lines[0])
    records = [line.split() for line in lines[2:2 + natoms]]
    labels = [record[0] for record in records]
    coordinates_angstrom = np.array(
        [[float(value) for value in record[1:4]] for record in records]
    )
    return labels, coordinates_angstrom.reshape(-1) * units.Angst


def summarize(title, result):
    """Print frequencies and the first sampled initial condition."""
    nonzero = result["omega_au"] > 0.0
    print(f"\n{title}")
    print("frequencies / cm^-1:", result["omega_au"][nonzero] / units.inv_cm2Ha)
    print("first coordinates / bohr:", result["ics"][0]["q"])
    print("first momenta / a.u.:", result["ics"][0]["p"])


def main():
    labels, q_eq = read_xyz_geometry(REFERENCE_DIR / "water_equilibrium.xyz")

    # Route 1: read frequencies and Cartesian displacements from mode_*.xyz.
    from_modes = prepare_wigner_from_modes(
        labels=labels,
        q_eq=q_eq,
        mode_start=1,
        mode_end=3,
        mode_file_pattern=str(REFERENCE_DIR / "mode_{}.xyz"),
        mass_map=MASSES_DA,
        temperature=300.0,
        ntraj=5,
        seed=12345,
    )
    summarize("Sampling from normal-mode files", from_modes)

    # Route 2: inspect modes obtained by mass-weighting and diagonalizing a
    # Cartesian Hessian, then sample directly from that Hessian.
    hessian = np.loadtxt(REFERENCE_DIR / "water_hessian_au.txt")
    atomic_masses = [MASSES_DA[label] for label in labels]
    omega, d_cart, d_p, mass_cart = build_modes_from_hessian(
        hessian, atomic_masses
    )
    print("\nHessian has", np.count_nonzero(omega), "vibrational modes")
    print("canonical transform error:", np.max(np.abs(d_cart.T @ d_p - np.eye(9))))
    print("Cartesian masses / electron mass:", mass_cart)

    from_hessian = generate_wigner_from_hessian(
        q_eq=q_eq,
        hessian=hessian,
        masses=atomic_masses,
        temperature=300.0,
        ntraj=5,
        seed=12345,
    )
    summarize("Sampling directly from the Hessian", from_hessian)


if __name__ == "__main__":
    main()

