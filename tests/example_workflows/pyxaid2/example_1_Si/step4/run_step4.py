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


params = {}

params["nstates"] = 4
params["nfiles"] = 30

# Normal example
params["data_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/pyxaid2/example_1_Si/step3/res_setup1/"]

# Pretend we have 2 directories:
#params["data_set_paths"].append("/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/pyxaid2/example_1_Si/step3/res_setup1/")

params["Hvib_re_prefix"] = "Hvib_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_"; params["Hvib_im_suffix"] = "_im"

# General simulaiton parameters
params["T"] = 300.0               # Temperature, K
params["ntraj"] = 10              # how many stochastic trajectories
params["sh_method"] = 1           # 0 - MSSH, 1 - FSSH
params["decoherence_method"] = 0  # 0 - no decoherence, 1 - decoherence (ID-A), 2 - MSDM, 3 - DISH
params["dt"] = 41.3413            # Nuclear dynamics integration timestep, in a.u.
params["nsteps"] = 10             # The length of the NA-MD trajectory
params["Boltz_opt"] = 3           # Option for the frustrated hops acceptance/rejection
params["istate"] = 3              # The index of the starting excited state (indexing from 0)
params["init_times"] = [0,15]     # starting points for sub-trajectories
params["outfile"] = "_out.txt"    # output file

step4.run(params)

