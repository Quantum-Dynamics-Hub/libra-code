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

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/math_random")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cygrandom import *
    from cygdyn import *
    from cyghamiltonian import *
    from cyglinalg import *


# Fisrt, we add the location of the library to test to the PYTHON path
#if sys.platform=="cygwin":
#    from cyglibra_core import *
#elif sys.platform=="linux" or sys.platform=="linux2":
#    from liblibra_core import *
#from libra_py import *





print "\nTest 2: Now we will be using model Hamiltonian"
# First, create the Hamiltonian
ham = Hamiltonian_Model(1)  # DAC
ham.set_rep(1)  # adiabatic


# Electronic - 2 levels, starting at 0-th state
el = Electronic(2,0)

# Nuclear DOFs
mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = -10.0
mol.p[0] = 20.0


f = open("_tsh.txt","w")
dt = 1.0

# Initialization
ham.set_v([ mol.p[0]/mol.mass[0] ])
#epot = compute_potential_energy(mol, el, ham, 1)  # 0 - MF forces, 1 - FSSH (state-resolved) forces
#compute_forces(mol, el, ham, 1) # 0 - MF, 1 - FSSH

epot = compute_forces(mol, el, ham, 1)  # 0 - MF forces, 1 - FSSH (state-resolved) forces


use_boltz_factor = 0
T = 300.0
rnd = Random()

for i in xrange(2500):

    # el-dyn
    el.propagate_electronic(0.5*dt, ham)

    # Nuclear dynamics
    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)

    ham.set_v([ mol.p[0]/mol.mass[0] ])
    epot = compute_forces(mol, el, ham, 1) 
    #epot = compute_potential_energy(mol, el, ham, 1)  # 0 - MF forces, 1 - FSSH
    #compute_forces(mol, el, ham, 1)

    mol.propagate_p(0.5*dt)

    # el-dyn
#    ham.set_v([ mol.p[0]/mol.mass[0] ])
    el.propagate_electronic(0.5*dt, ham)    
#    ham.set_v([ mol.p[0]/mol.mass[0] ])

#    ekin = compute_kinetic_energy(mol)


    # Now, incorporate surface hop
    g = MATRIX(2,2)

    # Just choose the TSH scheme below
#    compute_hopping_probabilities_mssh(mol, el, ham, g, dt, use_boltz_factor, T)
    compute_hopping_probabilities_fssh(mol, el, ham, g, dt, use_boltz_factor, T)
#    compute_hopping_probabilities_gfsh(mol, el, ham, g, dt, use_boltz_factor, T)

 
    do_rescaling = 1
    do_reverse = 1
    rep = 1 # adiabatic
    ksi = rnd.uniform(0.0, 1.0)

    #print i, ksi
    state_old = el.istate
    p_old = mol.p[0]
    f_old = mol.f[0]
    
    el.istate = hop(el.istate, mol, ham, ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary

    ekin = compute_kinetic_energy(mol)
    epot = compute_forces(mol, el, ham, 1)  
#    epot = compute_potential_energy(mol, el, ham, 1)


    if el.istate is not state_old:
        print "old p = ", p_old, "Ekin_old = ", 0.5 *p_old**2 / mol.mass[0]  , "Epot_old = ", ham.H(state_old,state_old).real, "etot_old = ", 0.5 *p_old**2 / mol.mass[0] + ham.H(state_old,state_old).real
        print "new p = ",mol.p[0], "Ekin_new = ", 0.5 *mol.p[0]**2 / mol.mass[0], "Epot_new = ", ham.H(el.istate,el.istate).real, "etot_new = ", 0.5 *mol.p[0]**2 / mol.mass[0] + ham.H(el.istate,el.istate).real
        print "ekin = ", ekin, " epot = ", epot
        print "f_old = ", f_old, "f_new = ", mol.f[0], " dHdR = ", ham.dHdq(el.istate,el.istate, 0)


    print i, ekin, epot, ekin+epot

    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f  |c0|^2= %8.5f  |c1|^2= %8.5f  Re|c01|= %8.5f istate= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot, el.rho(0,0).real, el.rho(1,1).real, el.rho(0,1).real, el.istate))

      
f.close()

