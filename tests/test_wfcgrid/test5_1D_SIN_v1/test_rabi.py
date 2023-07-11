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
wfc = Wfcgrid(-20.0, 30.0, 0.025, 2)

print "\nTest 3: Init wfc"
#   void init_wfc_1D(double x0, double px0, double dx, int init_state); // 1D
wfc.init_wfc_1D(0.0, 5.0, 1.0, 1)

#print "\nTest 4: Print snap"
wfc.print_wfc_1D("res/wfc",0, 0)
wfc.print_wfc_1D("res/wfc",0, 1)
wfc.print_reci_wfc_1D("res/reci_wfc",0, 0)
wfc.print_reci_wfc_1D("res/reci_wfc",0, 1)



# Initialize Hamiltonian: ECWR in adiabatic representation
x = 0.0
v = 1.00
ham = Hamiltonian_Model(7)
ham.set_rep(0)  # diabatic
E0 = -0.001
E1 = 0.001
V01= 0.001 # * math.sqrt(10.0)
D = 0.0019
L = 1.0
ham.set_params([E0, E1, V01, D, L ])

#ham.set_q([x])  
ham.set_v([v])  
#ham.compute()


print "\nTest 5: Update potential"
wfc.update_potential_1D(ham)

dt = 0.1
print "\nTest 6: Update propagator"
wfc.update_propagator_1D(0.5*dt, 1.0)  # this is important because we are using exp(-0.5*dt*H_loc)...


print "\nTest 7: Update propagator_K"
#wfc.update_propagator_K_1D(dt, 1.0)    # ... together with exp(-dt*H_non-loc)
wfc.update_propagator_K_1D(dt, 2000.0)    # ... together with exp(-dt*H_non-loc)


print "\nTest 8: Printing different info"
wfc.print_ham_1D("H_00.txt",0,0)
wfc.print_ham_1D("H_01.txt",0,1)
wfc.print_ham_1D("H_11.txt",1,1)

wfc.print_expH_1D("expH_00.txt",0,0)
wfc.print_expH_1D("expH_01.txt",0,1)
wfc.print_expH_1D("expH_11.txt",1,1)

wfc.print_expK_1D("expK_00.txt",0)
wfc.print_expK_1D("expK_11.txt",1)


#sys.exit(0)


print "\nTest 9: Propagate"
# Compute dynamics
for i in range(1,201):  # time steps
    print i
    for j in range(0,250):  # time steps
        wfc.propagate_exact_1D(0)

    wfc.print_wfc_1D("res/wfc", i, 0)
    wfc.print_wfc_1D("res/wfc", i, 1)
    wfc.print_reci_wfc_1D("res/reci_wfc", i, 0)
    wfc.print_reci_wfc_1D("res/reci_wfc", i, 1)
    wfc.print_populations_1D("populations.dat",i)


