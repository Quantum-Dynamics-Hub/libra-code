#*********************************************************************************
#* Copyright (C) 2018 Wei Li and Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************
#
# This script shows how to run Quasistochastic Hamiltonian NA-MD
#

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *
import libra_py.workflows.nbra.qsh as qsh


params = {}

params["dt"] = 41.0
params["nfreqs"] = 25 # number of frequency included 
params["dw"] = 1.0
params["wspan"] = 3000.0
params["filename"] = "_spectr_"
params["logname"] = "_log.txt"

params["nstates"] = 2
params["nfiles"] = 1000  # even though there are actually 1000 files
params["data_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/qsh/step2/res/"]
params["Hvib_re_prefix"] = "Hvib_"; params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "Hvib_"; params["Hvib_im_suffix"] = "_im"

params["do_output"] = True
params["nsteps"] = 500 # total length of quasi-stochastical Ham
params["output_set_paths"] = ["/mnt/c/cygwin/home/Alexey-user/Programming/Project_libra/tests/example_workflows/qsh/step3/res_qsh/"]
params["qsh_Hvib_re_prefix"] = "Hvib_"; params["qsh_Hvib_re_suffix"] = "_re"
params["qsh_Hvib_im_prefix"] = "Hvib_"; params["qsh_Hvib_im_suffix"] = "_im"


qsh.run(params)
