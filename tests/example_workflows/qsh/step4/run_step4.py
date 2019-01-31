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
import libra_py.workflows.nbra.qsh as qsh


#=========== QSH run (on-the-fly) =========================

qsh_params = {}
qsh_params["dt"] = 41.0
qsh_params["dw"] = 1.0
qsh_params["wspan"] = 3000.0

qsh_params["filename"] = "_spectr_"
qsh_params["logname"] = "_log.txt"
qsh_params["do_output"] = False
qsh_params["T"] = 300.0               # Temperature, K

qsh_params["nfreqs"] = 6 # number of frequency included 
qsh_params["norbitals"] = 56
qsh_params["active_space"] = [35,36]
qsh_params["nstates"] = 2
qsh_params["nfiles"] = 2000  # even though there are actually 1000 files
qsh_params["ntraj"] = 100
qsh_params["nsteps"] = 5000000 # total length of quasi-stochastical Ham
qsh_params["scale"] = 0.5             # scale factor for the Hvib matrix elements   
qsh_params["data_set_paths"] = ["/home/weili/QSH/res-mapbi/res/"]
qsh_params["Hvib_re_prefix"] = "0_Ham_"; qsh_params["Hvib_re_suffix"] = "_re"
qsh_params["Hvib_im_prefix"] = "0_Ham_"; qsh_params["Hvib_im_suffix"] = "_im"

qsh_params["sh_method"] = 1           # 0 - MSSH, 1 - FSSH
qsh_params["decoherence_method"] = 2  # 0 - no decoherence, 1 - decoherence (ID-A), 2 - MSDM, 3 - DISH
qsh_params["Boltz_opt"] = 3           # Option for the frustrated hops acceptance/rejection
qsh_params["istate"] = 1              # The index of the starting excited state (indexing from 0)
qsh_params["init_times"] = [0]  

qsh_params["on-the-fly-qsh"] = True
qsh_params["outfile"] = "_out_qsh_fly.txt"    # output file

qsh_params["set_decoh"] = True  # whether to set decoherence time manually
qsh_params["decoh_time"] = 5.6 # in fs


Hvib=step4.get_Hvib2(qsh_params)
step4.transform_data(Hvib,qsh_params)
step4.run(qsh.run(Hvib,qsh_params),qsh_params)


