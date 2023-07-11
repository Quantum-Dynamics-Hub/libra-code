#*********************************************************************************
#* Copyright (C) 2017-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/


import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

from common import compute_statistics, compute_statistics2, lattice, lattice_p, compute_rmsd
from Hamiltonian import compute_model
from md import run_nve, run_nvt



# ===== Define the system ==========
Nx, Ny, Nz = 5, 5, 5
nnucl = Nx * Ny * Nz
ndof, ntraj = 3*nnucl, 1

# Dynamical variables and system-specific properties
dx, dy, dz = 2.0, 2.0, 2.0
sx, sy, sz = 0.1, 0.1, 0.1
px, py, pz = 0.0, 0.0, 0.0
spx, spy, spz = 0.1, 0.1, 0.1

rnd = Random()
q = MATRIX(ndof,ntraj);  lattice(q, Nx,Ny,Nz, dx, dy, dz, sx, sy, sz, rnd)
p = MATRIX(ndof,ntraj);  lattice_p(p, Nx,Ny,Nz, px, py, pz, spx, spy, spz, rnd)

                         
iM = MATRIX(ndof,1);     
for i in xrange(ndof):
    iM.set(i,0, 1.0/2000.0)


# ===== Define the simulation parameters ==========
# Model parameters 
params = {}
params["x0"], params["k"], params["D"], params["V"], params["omega"] = 2.0, 0.05, -0.1, 0.05, 0.25

# Simulation and output parameters
params["nsteps"] = 1000
params["out_energy"] = "_output-anneal.txt"
params["out_phase_space"] = "_phase_space-anneal.txt"
params["out_positions"] = "_pos_space-anneal.txt"
params["label"] = ["H"]*nnucl
params["trajectory_index"] = 0
params["trajectory_filename"] = "_traj-anneal.xyz"

# Production run
# In case, we need a thermostat:
params["Temperature"] = 300.0
params["Q"] = 100.0
params["thermostat_type"] = "Nose-Hoover"
params["nu_therm"] = 0.001
params["NHC_size"] = 5


params["sigma"] = 2.0
params["epsilon"] = 0.01

# Annealing
params["dt"] = 10.0
params["nsteps"] = 50
for i in xrange(10):
    Q, P = run_nve(ndof, ntraj, q, p, iM, compute_model, params)

    # Reset velocities
    q = MATRIX(Q[-1])    
    lattice_p(p, Nx,Ny,Nz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rnd)


params["out_energy"] = "_output.txt"
params["out_phase_space"] = "_phase_space.txt"
params["out_positions"] = "_pos_space.txt"
params["trajectory_filename"] = "_traj.xyz"

params["dt"] = 10.0
params["nsteps"] = 2000
lattice_p(p, Nx,Ny,Nz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, rnd)
Q, P = run_nvt(ndof, ntraj, q, p, iM, compute_model, params)

idof = 0
minx = -1.0
maxx = 1.0
dx = 0.01

compute_statistics(Q, idof, minx, maxx, dx, "_density_q.txt")        
compute_statistics(P, idof, minx, maxx, dx, "_density_p.txt")        

compute_statistics2(Q, idof, 1, minx, maxx, dx, "_density_q-first.txt")        
compute_statistics2(P, idof, 1, minx, maxx, dx, "_density_p-first.txt")        


vel = []
for p in P:
    vel.append(p.col(0))
acf_matrix.recipe1(vel, params["dt"]*0.02419, 5000.0, 1.0, "_acf_vel.txt", "_spectrum_vel.txt", 0)


compute_rmsd(Q, params["dt"], "_rmsd.txt")

