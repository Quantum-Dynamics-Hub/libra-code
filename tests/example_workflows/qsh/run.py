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

from libra_py.workflows.qsh import *

params = {}
params["dt"] = 41.0
params["T"] = 300.0 

params["nfiles"] = 300 # number of direct Ham files we want to include 
params["nsteps"] = 20000 # total length of quasi-stochastical Ham
params["Nfreqs"] = 1 # number of nuclear modes included 
params["ntraj"] = 100 #number of stochastic SH 

# 35 - HOMO
# 36 - LUMO
params["norbitals"] = 56
params["active_space"] = [35,36]   # L -> H
params["istate"] = 1 # initial state

params["set_deco"] = 1 # whether to set decoherence time manually, 0 - no; 1 - yes
params["deco_time"]  =  6.5 # dephasing time obtained from the direct vibronic Ham files, in fs


params["rt"] = "res/" #directory containning the direct Ham files
params["Hvib_re_prefix"] = "0_Ham_"
params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "0_Ham_"
params["Hvib_im_suffix"] = "_im"

qstochastic_Ham.run(params)
