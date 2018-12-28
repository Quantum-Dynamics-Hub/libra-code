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

# General simulaiton parameters
params["nstates"] = 2
params["T"] = 300.0               # Temperature, K
params["ntraj"] = 100             # how many stochastic trajectories
params["sh_method"] = 1           # 0 - MSSH, 1 - FSSH
params["decoherence_method"] = 3  # 0 - no decoherence, 1 - decoherence (ID-A), 2 - MSDM, 3 - DISH
params["dt"] = 41.3413            # Nuclear dynamics integration timestep, in a.u.
params["nsteps"] = 300            # The length of the NA-MD trajectory
params["Boltz_opt"] = 3           # Option for the frustrated hops acceptance/rejection
params["istate"] = 1              # The index of the starting excited state (indexing from 0)
params["init_times"] = [0]        # starting points for sub-trajectories
params["Hvib_re_prefix"] = "Hvib_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_"; params["Hvib_im_suffix"] = "_im"


#=========== Normal run =========================
params["nfiles"] = 500
params["data_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/qsh/step2/res/"]
params["outfile"] = "_out.txt"    # output file

step4.run(params)


#=========== QSH run =========================
params["nfiles"] = 500
params["data_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/qsh/step3/res_qsh/"]
params["outfile"] = "_out_qsh.txt"    # output file

step4.run(params)


#=========== QSH run (on-the-fly) =========================

qsh_params = {}

qsh_params["dt"] = 41.0
qsh_params["nfreqs"] = 25 # number of frequency included 
qsh_params["dw"] = 1.0
qsh_params["wspan"] = 3000.0
qsh_params["filename"] = "_spectr_"
qsh_params["logname"] = "_log.txt"
qsh_params["nstates"] = 2
qsh_params["nfiles"] = 1000  # even though there are actually 1000 files
qsh_params["data_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/qsh/step2/res/"]
qsh_params["Hvib_re_prefix"] = "Hvib_"; params["Hvib_re_suffix"] = "_re"
qsh_params["Hvib_im_prefix"] = "Hvib_"; params["Hvib_im_suffix"] = "_im"
qsh_params["nsteps"] = 500 # total length of quasi-stochastical Ham

qsh_params["do_output"] = False

params["on-the-fly-qsh"] = True
params["qsh-params"] = qsh_params
params["outfile"] = "_out_qsh_fly.txt"    # output file

step4.run(params)


