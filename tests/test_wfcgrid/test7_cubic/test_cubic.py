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
nsteps = 20 # steps per snap
dt = 1.0  # step size a.u.
dx = 0.1
istate = 0

cases = [1]
x0 = [-0.2] #, -0.3, -0.4]
px0 = [0.0] #, 0.0, 0.0]



# Grid
wfc = Wfcgrid(-20.0, 20.0, 0.005, 1)

for A in cases:

    a = cases.index(A)

    # Hamiltonian
    ham = Hamiltonian_Model(8)
    ham.set_rep(0)  # diabatic
    ham.set_v([px0[a] / mass])  

    # Potential
    wfc.update_potential_1D(ham)
    wfc.print_ham_1D("H_00.txt",0,0)

    # Propagators
    wfc.update_propagator_1D(0.5*dt, mass)  # this is important because we are using exp(-0.5*dt*H_loc)...
    wfc.update_propagator_K_1D(dt, mass)    # ... together with exp(-dt*H_non-loc)



    os.system("mkdir res%i" % A)

    # Init wfc
    wfc.init_wfc_1D(x0[a], px0[a], dx, istate)


    # Propagate wfc
    f = open("abs.txt", "w")
    f.close()

    cum_l = 0.0
    cum_r = 0.0

    for i in range(0, nsnaps):  # time steps

        wfc.print_wfc_1D("res%i/wfc" % A, i, 0)
        wfc.print_reci_wfc_1D("res%i/reci_wfc" % A, i, 0)

        wfc.print_populations_1D("populations_%i.dat" % A,i)

        for j in range(0,nsteps):  # time steps
            wfc.propagate_exact_1D(0)

            pop_l, pop_r = wfc.absorb_1D(10.0)
            cum_l = cum_l + pop_l[0]
            cum_r = cum_r + pop_r[0]

        e_kin = wfc.e_kin_1D(mass)
        e_pot = wfc.e_pot_1D()

        f = open("abs.txt", "a")
        f.write("%10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n" % (dt*i*nsteps, cum_l, cum_r, e_kin, e_pot, e_kin+e_pot))
        f.close()





