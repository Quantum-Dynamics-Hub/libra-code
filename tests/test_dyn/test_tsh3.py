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

######################################################################
#
# Basic test for surface hopping schemes - now using object of Ensemble
# type. We can choose TSH methods: FSSH, MSSH, GFSH;  we can choose
# Tully's models and initial conditions. Computed populations are 
# resolved only by state, not like reflection/transmission
#
######################################################################


# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



print "\nTest 2: Set parameters"
ntraj = 150
use_boltz_factor = 0
T = 300.0
rnd = Random()
do_rescaling = 1 
do_reverse = 1
rep = 1 # adiabatic
ham_indx = 1
nel = 2
nnucl = 1
dt = 1.0


P_list = []
for i in xrange(60):
    P_list.append(1.0 + 0.5*i)

#P_list = [1.0, 10.0, 25.0, 50.0]


f = open("p_scatt.txt","w")
f.close()

for P in P_list:

    print "P = ", P
    print "Step 1: Initialize ensemble and Hamiltonians"
    ens = Ensemble(ntraj, nel, nnucl) 
    ens.ham_set_ham("model", ham_indx)
    ens.ham_set_rep(rep)

    # Initialization
    for i in xrange(ntraj):
        ens.mol[i].mass[0] = 2000.0
        ens.mol[i].q[0] = -10.0
        ens.mol[i].p[0] = P        
        ens.ham_set_v(i, [ens.mol[i].p[0]/ens.mol[i].mass[0]] )

    epot = compute_potential_energy(ens, 1)  # 0 - MF forces, 1 - FSSH
    compute_forces(ens, 1)

    #---------------------- Propagation -------------------------------
    print "Step 2: Propagation"
    for t in xrange(2500):
        ens.el_propagate_electronic(0.5*dt)

        # Nuclear dynamics
        ens.mol_propagate_p(0.5*dt)
        ens.mol_propagate_q(dt)
        ens.ham_set_v()

        epot = compute_potential_energy(ens, 1)  # 0 - MF forces, 1 - FSSH
        compute_forces(ens, 1)

        ens.mol_propagate_p(0.5*dt)

        # el-dyn
        ens.el_propagate_electronic(0.5*dt)    
        ens.ham_set_v()

        ekin = compute_kinetic_energy(ens)


        # Now, incorporate surface hop
        g = MATRIX(2,2)

        for i in xrange(ntraj):
        # Just choose the TSH scheme below
            #compute_hopping_probabilities_mssh(ens, i, g, dt, use_boltz_factor, T)
            compute_hopping_probabilities_fssh(ens, i, g, dt, use_boltz_factor, T)
            #compute_hopping_probabilities_gfsh(ens, i, g, dt, use_boltz_factor, T) 

            ksi = rnd.uniform(0.0, 1.0)
            ens.el[i].istate = hop(ens.el[i].istate, ens, i , ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary

    #---------------------- Propagation is done ------------------------------

    print "Step 3: Analysis and print"
    se_pops = ens.se_pop()
    sh_pops = ens.sh_pop()

    f = open("p_scatt.txt","a")
    f.write("P= %8.5f  se_pop(0)= %8.5f  se_pop(1)= %8.5f  sh_pop(0)= %8.5f sh_pop(1)= %8.5f\n" % (P, se_pops[0], se_pops[1], sh_pops[0], sh_pops[1] ))
    f.close()

print "time propagation is done"      
print "file is closed"


