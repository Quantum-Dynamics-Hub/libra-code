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

from common import sample, compute_statistics, compute_statistics2, chain_x, chain_px
from Hamiltonian import compute_model
from md import run_nve



# ===== Define the system ==========
nnucl = 5
ndof, ntraj = 3*nnucl, 5

# Dynamical variables and system-specific properties
dx = 2.1
sigma_q = 0.01
mean_p = 0.0
sigma_px = 0.1

rnd = Random()
q = MATRIX(ndof,ntraj);  chain_x(q, dx, sigma_q, 0.0, 0.0, rnd)
p = MATRIX(ndof,ntraj);  chain_px(p, mean_p, 0.0, 0.0, sigma_px, 0.0, 0.0,  rnd)
iM = MATRIX(ndof,1);     
for i in xrange(ndof):
    iM.set(i,0, 1.0/100.0)


# ===== Define the simulation parameters ==========
# Model parameters 
params = {}
params["x0"], params["k"], params["D"], params["V"], params["omega"] = 2.0, 0.05, -0.1, 0.05, 0.25

# Simulation and output parameters
params["dt"] = 10.0
params["nsteps"] = 2000
params["out_energy"] = "_output.txt"
params["out_phase_space"] = "_phase_space.txt"
params["out_positions"] = "_pos_space.txt"
params["label"] = ["H"]*nnucl
params["trajectory_index"] = 0
params["trajectory_filename"] = "_traj.xyz"


Q, P = run_nve(ndof, ntraj, q, p, iM, compute_model, params)

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
