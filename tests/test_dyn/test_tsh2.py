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


print "\nTest 2: Initialize ensemble"
ntraj = 50
#sys.exit(0)
ens = Ensemble(ntraj, 2, 1)  # num_of_traj,  num_el, num_nucl


print "\nTest 3: Now we will be using model Hamiltonian"

# One way to setup Hamiltonian
ens.ham_set_ham("model", 1)
ens.ham_set_rep(1)

for i in xrange(ntraj):
    print i
# Alternative way:
#    ens.ham_set_ham(i, "model", 1)
#    ens.ham_set_rep(i, 1)  # adiabatic

    # Nuclear DOFs
    ens.mol[i].mass[0] = 2000.0
    ens.mol[i].q[0] = -10.0
    ens.mol[i].p[0] = 20.0


print "hy-hy"
#sys.exit(0)




f = open("_tsh.txt","w")
dt = 1.0


# Initialization
#ens.ham_set_v()

for i in xrange(ntraj):
    print i
    ens.ham_set_v(i, [ens.mol[i].p[0]/ens.mol[i].mass[0]] ) 

#    epot = compute_potential_energy(ens.mol[i], ens.el[i], ens.ham[i], 1)  # 0 - MF forces, 1 - FSSH (state-resolved) forces
#    compute_forces(ens.mol[i], ens.el[i], ens.ham[i], 1) # 0 - MF, 1 - FSSH

print "hy-hy-hy"


use_boltz_factor = 0
T = 300.0
rnd = Random()
do_rescaling = 1 
do_reverse = 1
rep = 1 # adiabatic


print "Generated random"
#sys.exit(0)


for t in xrange(2500):

    print "t = ",t
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
        compute_hopping_probabilities_mssh(ens, i, g, dt, use_boltz_factor, T)
        #compute_hopping_probabilities_fssh(ens, i, g, dt, use_boltz_factor, T)
        #compute_hopping_probabilities_gfsh(ens, i, g, dt, use_boltz_factor, T) 

        ksi = rnd.uniform(0.0, 1.0)
        ens.el[i].istate = hop(ens.el[i].istate, ens, i , ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary


    se_pops = ens.se_pop()
    sh_pops = ens.sh_pop()

    f.write("t= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  se_pop(0)= %8.5f  se_pop(1)= %8.5f  sh_pop(0)= %8.5f sh_pop(1)= %8.5f\n" % (t, ekin, epot, ekin+epot, se_pops[0], se_pops[1], sh_pops[0], sh_pops[1] ))

print "time propagation is done"      

f.close()
print "file is closed"


