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
# Wfcgrid::Wfcgrid(double minx_, double maxx_, double dx_, double miny_, double maxy_, double dy_, int nstates_)
wfc = Wfcgrid(-10.0, 15.0, 0.1, -10.0, 15.0, 0.1, 2)

print "\nTest 3: Init wfc"
# void Wfcgrid::init_wfc_2D(double x0, double y0, double px0, double py0, double dx0, double dy0, int init_state)
wfc.init_wfc_2D(0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 1)

#print "\nTest 4: Print snap"
wfc.print_wfc_2D("wfc",0, 0)


# Initialize Hamiltonian: SIN in diabatic representation
x = 0.0
v = 1.0
ham = Hamiltonian_Model(200)
ham.set_rep(0)  # diabatic

E0 = -0.001
E1 = 0.001
V01= 0.001 # * math.sqrt(10.0)
D = 0.0019
Lx = 1.0
Ly = 1.0


ham.set_params([E0, E1, V01, D, Lx, Ly ])
ham.set_q([x])  
ham.set_v([v])  
ham.compute()


print "\nTest 5: Update potential"
wfc.update_potential_2D(ham)

dt = 0.1
print "\nTest 6: Update propagator"
wfc.update_propagator_2D(0.5*dt, 2000.0)  # this is important because we are using exp(-0.5*dt*H_loc)...

print "\nTest 7: Update propagator_K"
wfc.update_propagator_K_2D(dt, 2000.0)    # ... together with exp(-dt*H_non-loc)


print "\nTest 8: Propagate"
# Compute dynamics
for i in range(0,201):  # time steps
    print i
    for j in range(0,250):  # time steps
        wfc.propagate_exact_2D(2)
    wfc.print_wfc_2D("res/wfc", i, 0)
    wfc.print_wfc_2D("res/wfc", i, 1)
    wfc.print_populations_2D("populations.dat",i)


