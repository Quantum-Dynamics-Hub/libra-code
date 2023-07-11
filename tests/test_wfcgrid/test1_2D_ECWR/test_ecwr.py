#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



print "\nTest 2: Init grid"
wfc = Wfcgrid(-20.0, 20.0, 0.1, -10.0, 10.0, 1.0, 2)

print "\nTest 3: Init wfc"
wfc.init_wfc_2D(-10.0, 0.0, 20.0, 0.0, 0.1, 0.01, 0)

#print "\nTest 4: Print snap"
wfc.print_wfc_2D("wfc",0, 0)


# Initialize Hamiltonian: ECWR in adiabatic representation
x = -10.0
v = 1.0
ham = Hamiltonian_Model(2)
ham.set_rep(0)  # diabatic
ham.set_q([x])  
ham.set_v([v])  
ham.compute()


print "\nTest 5: Update potential"
wfc.update_potential_2D(ham)

dt = 1.0
print "\nTest 6: Update propagator"
wfc.update_propagator_2D(0.5*dt, 2000.0)  # this is important because we are using exp(-0.5*dt*H_loc)...

print "\nTest 7: Update propagator_K"
wfc.update_propagator_K_2D(dt, 2000.0)    # ... together with exp(-dt*H_non-loc)


print "\nTest 8: Propagate"
# Compute dynamics
for i in range(0,25):  # time steps
    print i
    for j in range(0,50):  # time steps
        wfc.propagate_exact_2D(2)
    wfc.print_wfc_2D("res/wfc", i, 0)
    wfc.print_wfc_2D("res/wfc", i, 1)




