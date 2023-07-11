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

#####################################################################################
#
# Here, we will show the usage of Hamiltonian_Extern class: Ehrenfest dynamics - this is 
# just the same as in the ../test_dyn folder
#
#####################################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



print "\nTest 2: Now we will be using model Hamiltonian"
# First, create the Hamiltonian
ham = Hamiltonian_Model(1)  # DAC - this is the actual Hamiltonian - it will mimic an external procedure to compute Hamiltonian
ham.set_rep(1)  # adiabatic

# Here goes the external Hamiltonian - the one which will actually be used
ham_ex = Hamiltonian_Extern(2,1)  # 2 - # of electronic DOF, 1 - # of nuclear DOF
ham_ex.set_rep(1)  # adiabatic
ham_ex.set_adiabatic_opt(0)  # use the externally-computed adiabatic electronic Hamiltonian and derivatives
ham_ex.set_vibronic_opt(0)  # use the externally-computed vibronic Hamiltonian and derivatives

# Actual matrices
ham_adi = MATRIX(2,2);  ham_ex.bind_ham_adi(ham_adi);
d1ham_adi = MATRIXList()

tmp = MATRIX(2,2)
d1ham_adi.append(tmp);  ham_ex.bind_d1ham_adi(d1ham_adi);

ham_vib = CMATRIX(2,2);  ham_ex.bind_ham_vib(ham_vib);




# Electronic - 2 levels, starting at 0-th state
el = Electronic(2,0)

# Nuclear DOFs
mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = -10.0
mol.p[0] = 20.0


f = open("_ehrenfest_extern.txt","w")
dt = 1.0


# Update matrices
ham.set_v([ mol.p[0]/mol.mass[0] ])
ham.set_q([ mol.q[0] ])
ham.compute()
for i in [0,1]:
    for j in [0,1]:
        ham_adi.set(i,j, ham.H(i,j).real)
        ham_vib.set(i,j, ham.Hvib(i,j))
        d1ham_adi[0].set(i,j, ham.dHdq(i,j,0).real)


# Initialization
ham_ex.set_v([ mol.p[0]/mol.mass[0] ])
epot = compute_potential_energy(mol, el, ham_ex, 0)  # 0 - MF forces
compute_forces(mol, el, ham_ex, 0)


#sys.exit(0)

for t in xrange(5000):

    # el-dyn
    el.propagate_electronic(0.5*dt, ham_ex)

    # Nuclear dynamics
    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)



    # Update matrices
    ham.set_v([ mol.p[0]/mol.mass[0] ])
    ham.set_q([ mol.q[0] ])
    ham.compute()
    for i in [0,1]:
        for j in [0,1]:
            ham_adi.set(i,j, ham.H(i,j).real)
            ham_vib.set(i,j, ham.Hvib(i,j))
            d1ham_adi[0].set(i,j, ham.dHdq(i,j,0).real)


    ham.set_v([ mol.p[0]/mol.mass[0] ])
    epot = compute_potential_energy(mol, el, ham_ex, 0)  # 0 - MF forces
    compute_forces(mol, el, ham_ex, 0)

    mol.propagate_p(0.5*dt)

    # el-dyn
    el.propagate_electronic(0.5*dt, ham_ex)

    ekin = compute_kinetic_energy(mol)

    f.write("t= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f  |c0|^2= %8.5f  |c1|^2= %8.5f  Re|c01|= %8.5f\n" % (t, mol.q[0], mol.p[0], ekin, epot, ekin+epot, el.rho(0,0).real, el.rho(1,1).real, el.rho(0,1).real))

      
f.close()

