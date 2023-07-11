#*********************************************************************************
#* Copyright (C) 2015-2016 Brendan A. Smith, Alexey V. Akimov
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


# Here we initialize the grid and wavefunction
wfc = Wfcgrid(-5.0, 5.0, 0.01, 1)
wfc.init_wfc_1D(-1.0, 0.0, 0.1, 0)

# Here we make our Hamiltonian. 9 == symmetric double_well potential
ham = Hamiltonian_Model(9)
ham.set_rep(0)  # diabatic
ham.set_params([1.0])  ## Here, we set A = 1.0 for  V(q) = A( 0.25*x^4 - 0.5*x^2 )
a = open("dw_dia.txt","w")
for i in range(-500,500):
    q = 0.01*i
    ham.set_q([q])
    ham.set_v([0.0])
    ham.compute()
    a.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (q, ham.H(0,0).real, ham.H(1,1).real, ham.H(0,1).real ) )
    

ham.set_q([-1.0])
ham.set_v([0.0])
ham.compute()
wfc.update_potential_1D(ham)
dt = 1.0
wfc.update_propagator_1D(0.5*dt, 1835.0)  # this is important because we are using exp(-0.5*dt*H_loc)...
wfc.update_propagator_K_1D(dt, 1835.0)    # ... together with exp(-dt*H_non-loc)
                                          # 1835 = mass proton (a.u)
# Compute dynamics
for i in range(0,10):  # time steps
    for j in range(0,10):  # time steps
        wfc.propagate_exact_1D(0)
    wfc.print_wfc_1D("res/wfc", i, 0)   # Be sure to make a "res" directory


