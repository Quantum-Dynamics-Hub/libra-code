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
mol.p[0] = 10.0


f = open("_nucl_2.1.txt","w")
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

