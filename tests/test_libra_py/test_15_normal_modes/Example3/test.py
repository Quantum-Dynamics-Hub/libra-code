#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
###################################################################
#
# This example demonstrates how to obtain normal modes from LAMMPS
# 
###################################################################

import os
import math
import sys
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


import lammps 
lmp = lammps.lammps()

# Do the Hessian/dynmat calculations
D, H, R, M, typ  = LAMMPS_methods.compute_dynmat(lmp, "in.lammps", [1], 1e-3, 0)

print "Atom types are:", typ
print "Masses are:"; M.show_matrix()
print "Coordinates are:"; R.show_matrix()
#print "Hessian matrix is:"; H.show_matrix()
#print "Dynamical matrix is:"; D.show_matrix()

# Map atom types to atom names
E = []
PT = {1:"Au"}
for t in typ:
    E.append(PT[t])

# Visualizaiton and frequencies
params = {"nsteps":100, "nperiods":1, "scale":250.0, "print_modes":[5], "prefix":"_" }
normal_modes.compute_dynmat(R, D, M, E, params)



