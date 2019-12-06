#***********************************************************
# * Copyright (C) 2017-2018 Brendan A. Smith, Wei Li, and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys
import time
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.lz as lz

from libra_py import units
import libra_py.workflows.nbra.step4 as step4
import libra_py.workflows.nbra.lz as lz
import libra_py.workflows.nbra.decoherence_times as decoherence_times
from libra_py import data_conv
from libra_py import fit
from libra_py import influence_spectrum as infsp

params = {}

trajs = range(1)
#params["data_set_paths"] = []
#for itraj in trajs:
#    params["data_set_paths"].append( "/budgetdata/academic/alexeyak/brendan/Si_QD/H_capped/08nm/step2/traj"+str(itraj)+"/res/" )
#    params["data_set_paths"].append( "/budgetdata/academic/alexeyak/brendan/Si_QD/H_capped/08nm/step2/traj"+str(itraj)+"/res/" )

#case = 1   # 0.8 nm Si NC
case = 2   # 1.5 nm Si NC
#case = 3   # 2.2 nm Si NC

if case == 1:
    params["data_set_paths"] = ["res_0.8nm/"]
elif case == 2:
    params["data_set_paths"] = ["res_1.5nm/"]
elif case == 3:
    params["data_set_paths"] = ["res_2.2nm/"]




params["nfiles"]  = 100  # Ex) # of Hvib files to read for a given traj
params["Hvib_re_prefix"] = "hvib_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "hvib_"; params["Hvib_im_suffix"] = "_im"

# Set number of basis states
if case == 1:
    params["nstates"] = 30
elif case == 2:
    params["nstates"] = 100
elif case == 3:
    params["nstates"] = 300


# General simulaiton parameters
params["T"]                  = 300.0                 # Temperature, K
params["ntraj"]              = 100                   # how many stochastic trajectories
params["dt"]                 = 1.0*units.fs2au       # Nuclear dynamics integration timestep, in a.u.
params["nsteps"]             = params["nfiles"]      # The length of the NA-MD trajectory
params["istate"]             = 24                    # The index of the starting excited state (indexing from 0)
params["init_times"]         = [0]                   # starting points for sub-trajectories
params["do_output"]          = True                  # request to print the results into a file
params["do_return"]          = False                 # request to not store the date in the operating memory
params["gap_min_exception"]  = 0                     # set the minimal gap to zero if the extrapolated gap is negative
params["target_space"]       = 1                     # hop to all states


# For running NA-MD
Hvib = step4.get_Hvib2(params)      # get the Hvib for all data sets, Hvib is a lists of lists

init_time = params["init_times"][0]

# Compute average decoherence time over entire trajectory
tau, rates = decoherence_times.decoherence_times_ave(Hvib, [init_time], params["nfiles"]-init_time, 0)
avg_deco = tau/units.fs2au
#avg_deco.show_matrix()


#====================== One case =====================
# Looking on the "SH" populations - NBRA-TSH approach
params["Boltz_opt"]          = 1                     # Option for the frustrated hops acceptance/rejection
params["Boltz_opt_BL"]       = 0                     # Option to incorporate hte frustrated hops into BL probabilities
params["outfile"]            = "_out_TSH_.txt"       # output file
params["evolve_Markov"]      = False                 # don't propagate Markov
params["evolve_TSH"]         = True                  # Rely on propagating trajectories

start = time.time()
res = lz.run(Hvib, params)
end = time.time()
print("Time to run = ", end - start)



#====================== Another case =====================
# Looking on the "SE" populations - Markov chain approach
params["Boltz_opt"]          = 0                     # Option for the frustrated hops acceptance/rejection
params["Boltz_opt_BL"]       = 1                     # Option to incorporate hte frustrated hops into BL probabilities
params["outfile"]            = "_out_Markov_.txt"    # output file
params["evolve_Markov"]      = True                  # Rely on the Markov approach
params["evolve_TSH"]         = False                 # don't care about TSH

start = time.time()
res = lz.run(Hvib, params)
end = time.time()
print("Time to run = ", end - start)

