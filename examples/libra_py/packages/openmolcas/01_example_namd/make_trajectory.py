import os, sys
import numpy as np
import h5py
from liblibra_core import MATRIX, CMATRIX, CMATRIXList, Py2Cpp_int, Cpp2Py, Random
from libra_py import data_conv
from libra_py import units


#labels = ['C', 'C', 'C', 'H', 'H', 'H', 'H', 'O']
#labels = ['O', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
labels = ['Al', 'Al', 'Al']


F = h5py.File('FSSH2_/mem_data.hdf')
q = np.array(F['q/data'])
F.close()


nsteps, ntraj, ndof = q.shape
nat = ndof // 3


os.makedirs("trajectories", exist_ok=True)

for itraj in range(ntraj):
    filename = f"trajectories/traj_{itraj}.xyz"

    with open(filename, "w") as f:
        for istep in range(nsteps):
            f.write(f"{nat}\n")
            f.write(f"Trajectory {itraj}, step {istep}\n")

            for iat in range(nat):
                label = labels[iat]
                x = q[istep, itraj, 3 * iat + 0] / units.Angst
                y = q[istep, itraj, 3 * iat + 1] / units.Angst
                z = q[istep, itraj, 3 * iat + 2] / units.Angst

                f.write(
                    f"{label:<4s}"
                    f"{x:16.10f}"
                    f"{y:16.10f}"
                    f"{z:16.10f}\n"
                )


"""
os.system("mkdir trajectories")

for i in range(ntraj):
    #os.system(F"mkdir trajectories/traj_{i}.xyz")
    f = open(F"trajectories/traj_{i}.xyz", "w")
    for istep in range(nsteps):
        f.write(F"{nat}\n\n")
        for iat in range(nat):
            L = labels[iat]
            x = q[istep, i, 3*iat + 0] / units.Angst
            y = q[istep, i, 3*iat + 1] / units.Angst
            z = q[istep, i, 3*iat + 2] / units.Angst
            f.write(F"{L} {x}  {y}  {z}\n")
    f.close()
"""
