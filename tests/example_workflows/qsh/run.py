#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
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
params["nfiles"] = 2000 # number of direct direct Ham files we want to include 
params["nsteps"] = 20000000 # total length of stochastical direct Ham
params["T"] = 300.0 
params["Nfreqs"] = 6 # number of nuclear modes included in the calculation
params["ntraj"] = 100 

# 35 - HOMO
# 36 - LUMO
# 37 - LUMO+1
#params["active_space"] = [36,37]   # L +1 -> L
params["active_space"] = [35,36]   # L -> H
params["istate"] = 1 # initial state
params["rt"] = "/home/eric/PROJECT/Alexey3/res/" #directory containning the direct direct Ham files
params["deco_time"]  =  6.5 # dephasing time obtained from the direct vibronic Ham files, in fs

qstochastic_Ham.run(params)
