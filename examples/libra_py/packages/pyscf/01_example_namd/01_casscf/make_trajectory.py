"""Convert the saved Libra nuclear coordinates to XYZ trajectories."""

from pathlib import Path

import h5py
import numpy as np

from libra_py import units


ATOM_LABELS = ("Al", "Al", "Al")
INPUT_FILE = Path("FSSH2_/mem_data.hdf")
OUTPUT_DIRECTORY = Path("trajectories")


with h5py.File(INPUT_FILE, "r") as data_file:
    coordinates = np.asarray(data_file["q/data"])

nsteps, ntraj, ndof = coordinates.shape
natoms = ndof // 3
if ndof % 3 != 0 or len(ATOM_LABELS) != natoms:
    raise ValueError(
        f"Coordinate data contains {ndof} degrees of freedom, but "
        f"{len(ATOM_LABELS)} atom labels were provided"
    )

OUTPUT_DIRECTORY.mkdir(exist_ok=True)
for trajectory in range(ntraj):
    output_file = OUTPUT_DIRECTORY / f"traj_{trajectory}.xyz"
    with output_file.open("w", encoding="utf-8") as stream:
        for step in range(nsteps):
            stream.write(f"{natoms}\n")
            stream.write(f"Trajectory {trajectory}, step {step}\n")
            for atom, label in enumerate(ATOM_LABELS):
                xyz_angstrom = coordinates[step, trajectory, 3 * atom : 3 * atom + 3]
                xyz_angstrom = xyz_angstrom / units.Angst
                stream.write(
                    f"{label:<4s}"
                    f"{xyz_angstrom[0]:16.10f}"
                    f"{xyz_angstrom[1]:16.10f}"
                    f"{xyz_angstrom[2]:16.10f}\n"
                )
