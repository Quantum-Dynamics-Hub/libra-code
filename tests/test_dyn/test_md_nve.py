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
# This example demonstrates MD in NVE ensemble - for bigger system
#
######################################################################

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


#======================== Set force field and control parameters =============
ctx = Context()
ctx.set_path("md_nve")
param1 = 1
ctx.add("param1", param1)



print "\nTest 2: Now we will be using model Hamiltonian"
# First, create the Hamiltonian
ham = Hamiltonian_Atomistic()
ham.set_rep(1)  # adiabatic

ham.set_params([ff]) # !!!


# Electronic - 2 levels, starting at 0-th state
el = Electronic(2,0)

# Nuclear DOFs
N = 4 * (2*L)**3    # number of atoms
mol = Nuclear(3*N)  # number of DOFs

a = 5.256 * 1.889725989 # Angstrom to atomic units
# build fcc grid
I = 0
for i in range(-L,L):
    for j in range(-L,L):
        for k in range(-L,L):
            # Basis atom #1
            mol.q[I+0] = i*a  # x
            mol.q[I+1] = j*a  # y
            mol.q[I+2] = k*a  # z
            I = I + 3

            # Basis atom #2
            mol.q[I+0] = (i+0.5)*a  # x
            mol.q[I+1] = j*a  # y
            mol.q[I+2] = k*a  # z
            I = I + 3

            # Basis atom #3
            mol.q[I+0] = i*a  # x
            mol.q[I+1] = (j+0.5)*a  # y
            mol.q[I+2] = k*a  # z
            I = I + 3

            # Basis atom #4
            mol.q[I+0] = i*a  # x
            mol.q[I+1] = j*a  # y
            mol.q[I+2] = (k+0.5)*a  # z
            I = I + 3

for i in xrange(N):
    mol.mass[i] = 1821.5 * 39.95  # Ar amu to atomic units (m_e = 1)
    mol.p[i] = 0.0  # momenta


f = open("_md_nve.txt","w")
dt = 1.0

for i in xrange(5000):

    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)

    epot = compute_potential_energy(mol, el, ham, 1)  # 1 - FSSH forces
    compute_forces(mol, el, ham, 1)

    mol.propagate_p(0.5*dt)

    ekin = compute_kinetic_energy(mol)

    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot))

      
f.close()

