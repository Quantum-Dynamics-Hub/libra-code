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
params["nsteps"] = 200 # total length of quasi-stochastical Ham
params["Nfreqs"] = 1 # number of frequency included 
params["ntraj"] = 2000 #number of stochastic SH 

params["norbitals"] = 62
params["active_space"] = [30,31]   # L -> H
params["istate"] = 1 # initial state, start from zero

params["set_decoherence"] = -1 # how to set decoherence time , 
                               # -1 - don't include decoherence effects
                               # 0 - set decoherence time use Libra module; 
			       # 1 - set decoherence time manually

params["deco_time"]  =  6.5 # decoherence time obtained from the direct vibronic Ham files, in fs
                            # works only when set_decoherence = 1

params["time_inteval"] = 1000 # the QSH energies and couplings are printed out every time_interval steps
                              # QSH energies and couplings files will be stored in out directory

# input data
params["Hvib_re_prefix"] = "/home/weili/bzs/res/0_Ham_" # need to make res directory before running the script
params["Hvib_re_suffix"] = "_re"
params["Hvib_im_prefix"] = "/home/weili/bzs/res/0_Ham_"
params["Hvib_im_suffix"] = "_im"

# output data
params["qsh_Ham_prefix"] = "out/qsh_Ham_" # need to make out directory before running the script


qstochastic_Ham.run_namd(params)
