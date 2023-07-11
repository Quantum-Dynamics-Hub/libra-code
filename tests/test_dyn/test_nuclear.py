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


print "\nTest 2: Coupling nuclear dynamics to external forces - Harmonic oscillator"
def pot(q):
    k = 10.0
    return -k*q, 0.5*k*q*q

mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = 1.0
mol.p[0] = 0.0

f = open("_nucl.txt","w")

print "address mol = ", mol

dt = 1.0
for i in xrange(500):
    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)
    mol.f[0], epot = pot(mol.q[0])
    mol.propagate_p(0.5*dt)

    ekin = compute_kinetic_energy(mol)

    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot))

      
f.close()


print "\nTest 3: Coupling nuclear dynamics to external forces - Morse potential"
def pot2(q):
    k = 10.0
    beta = 1.1
    x = math.exp(-beta*q)
    return -k*(x-1.0)*(-beta*x), 0.5*k*(x*x - 2.0*x)

mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = 1.0
mol.p[0] = 0.0

f = open("_nucl2.txt","w")

print "address mol = ", mol

dt = 1.0
for i in xrange(2500):
    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)
    mol.f[0], epot = pot2(mol.q[0])
    mol.propagate_p(0.5*dt)

    ekin = compute_kinetic_energy(mol)

    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot))

      
f.close()



