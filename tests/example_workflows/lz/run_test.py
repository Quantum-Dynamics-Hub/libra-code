#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************
#
# This script shows how to run LZ NA-MD
#

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import libra_py.workflows.lz.lz as lz
import libra_py.units as units


params = {}
#===== Data description =====
params["norbitals"] = 2                # how many lines/columns in the file
params["active_space"] = [0, 1]        # which orbitals we care about (indexing starts with 0)
params["Hvib_re_prefix"] = "res/Hvib_" # prefixes of the files with real part of the Hamiltonian
params["Hvib_im_prefix"] = "res/Hvib_" # prefixes of the files with imaginary part of the Hamiltonian
params["Hvib_re_suffix"] = "_re"       # suffixes of the files with real part of the Hamiltonian
params["Hvib_im_suffix"] = "_im"       # suffixes of the files with imaginary part of the Hamiltonian
params["nsteps"] = 1000                # how many files to read

#===== Modeling params ===== 
params["dt"] = 41.0                    # [a.u.] time distance between the adjacent data points
params["ntraj"] = 1000                 # how many stochastic trajectories to use in the ensemble
params["istate"] = 1                   # index of the starting state (within those used in the active_space - see above)
params["outfile"] = "_test.txt"        # the name of the file, where all the results will be printed out
params["T"] = 600.0                    # [K] temperature of the simulation


lz.run_LZ(params)


