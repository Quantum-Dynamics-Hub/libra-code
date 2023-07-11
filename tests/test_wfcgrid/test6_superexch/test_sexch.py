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



# Parameters
mass = 2000.0
nsnaps = 200
nsteps = 1000 # steps per snap
dt = 0.1  # step size a.u.
x0 = -5.0
dx = 1.0
istate = 0

cases = range(0,50)
px0 = []
for i in cases:
    px0.append(0.0 + i*0.5)



# Grid
wfc = Wfcgrid(-20.0, 30.0, 0.025, 3)

for a in cases:

    # Hamiltonian
    ham = Hamiltonian_Model(4)
    ham.set_rep(0)  # diabatic
    ham.set_v([px0[a] / mass])  

    # Potential
    wfc.update_potential_1D(ham)
    wfc.print_ham_1D("H_00.txt",0,0)
    wfc.print_ham_1D("H_01.txt",0,1)
    wfc.print_ham_1D("H_02.txt",0,2)
    wfc.print_ham_1D("H_11.txt",1,1)
    wfc.print_ham_1D("H_12.txt",1,2)
    wfc.print_ham_1D("H_22.txt",2,2)

    # Propagators
    wfc.update_propagator_1D(0.5*dt, mass)  # this is important because we are using exp(-0.5*dt*H_loc)...
    wfc.update_propagator_K_1D(dt, mass)    # ... together with exp(-dt*H_non-loc)



    os.system("mkdir res%i" % a)

    # Init wfc
    wfc.init_wfc_1D(x0, px0[a], dx, istate)


    # Propagate wfc
    for i in range(0, nsnaps):  # time steps

        wfc.print_wfc_1D("res%i/wfc" % a, i, 0)
        wfc.print_wfc_1D("res%i/wfc" % a, i, 1)
        wfc.print_wfc_1D("res%i/wfc" % a, i, 2)
        wfc.print_reci_wfc_1D("res%i/reci_wfc" % a, i, 0)
        wfc.print_reci_wfc_1D("res%i/reci_wfc" % a, i, 1)
        wfc.print_reci_wfc_1D("res%i/reci_wfc" % a, i, 2)

        wfc.print_populations_1D("populations_%i.dat" % a,i)


        for j in range(0,nsteps):  # time steps
            wfc.propagate_exact_1D(0)




