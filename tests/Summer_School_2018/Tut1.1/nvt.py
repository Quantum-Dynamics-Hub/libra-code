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
from common import sample, compute_statistics, compute_statistics2
from Hamiltonian import compute_model
from md import run_nvt



# ===== Define the system ==========
nnucl, ntraj = 1, 25

# Dynamical variables and system-specific properties
mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, 0.1)
sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, 0.05)
mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, 0.0)
sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, 0.01)

rnd = Random()
q = MATRIX(nnucl,ntraj);  sample(q, mean_q, sigma_q, rnd)
p = MATRIX(nnucl,ntraj);  sample(p, mean_p, sigma_p, rnd)
iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/100.0)

# ===== Define the simulation parameters ==========
# Model parameters 
params = {}
params["x0"], params["k"], params["D"], params["V"], params["omega"] = 1.0, 0.1, -0.1, 0.05, 0.25

# Simulation and output parameters
params["dt"] = 10.0
params["nsteps"] = 2500
params["out_energy"] = "_output.txt"
params["out_phase_space"] = "_phase_space.txt"
params["out_positions"] = "_pos_space.txt"

# In case, we need a thermostat:
params["Temperature"] = 300.0
params["Q"] = 100.0
params["thermostat_type"] = "Nose-Hoover"
params["nu_therm"] = 0.01
params["NHC_size"] = 10


Q, P = run_nvt(nnucl, ntraj, q, p, iM, compute_model, params)

        
idof = 0
minx = -1.0
maxx = 1.0
dx = 0.01

compute_statistics(Q, idof, minx, maxx, dx, "_density_q.txt")        
compute_statistics(P, idof, minx, maxx, dx, "_density_p.txt")        

compute_statistics2(Q, idof, 1, minx, maxx, dx, "_density_q-first.txt")        
compute_statistics2(P, idof, 1, minx, maxx, dx, "_density_p-first.txt")        


