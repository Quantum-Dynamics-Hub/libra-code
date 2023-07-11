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



print "\nTest 2: Constructor #1: Default"
el = Electronic()
print "nstates = ", el.nstates
print "istate = ", el.istate
for i in range(0,el.nstates):
    print "q[%3i]= %8.5f , p[%3i]= %8.5f" % (i, el.q[i], i, el.p[i])


print "\nTest 3: Constructor #2: Arbitrary # of states, initial state is 0-th"
el = Electronic(3)
print "nstates = ", el.nstates
print "istate = ", el.istate
for i in range(0,el.nstates):
    print "q[%3i]= %8.5f , p[%3i]= %8.5f" % (i, el.q[i], i, el.p[i])


print "\nTest 4: Constructor #3: Arbitrary # of states, initial state is selected manually"
el = Electronic(3,2)
print "nstates = ", el.nstates
print "istate = ", el.istate
for i in range(0,el.nstates):
    print "q[%3i]= %8.5f , p[%3i]= %8.5f" % (i, el.q[i], i, el.p[i])


print "\nTest 5: Constructor #4: Arbitrary # of states, initial state is selected manually, also set the phase of the initial state"
el = Electronic(3,2, 0.45)
print "nstates = ", el.nstates
print "istate = ", el.istate
for i in range(0,el.nstates):
    print "q[%3i]= %8.5f , p[%3i]= %8.5f" % (i, el.q[i], i, el.p[i])



print "\nTest 6: Compute populations and densty matrices"
el = Electronic(3,1)
print "nstates = ", el.nstates
print "istate = ", el.istate
for i in range(0,el.nstates):
    print i, el.c(i)
    print "q[%3i]= %8.5f , p[%3i]= %8.5f, |c[%i]|^2= %8.5f" % (i, el.q[i], i, el.p[i], i, el.rho(i,i).real)

print "Density matrix:"
for i in range(0,el.nstates):
    tmp = ""
    for j in range(0,el.nstates):
        tmp = tmp + str(" %8.5f " % el.rho(i,j).real)
    print tmp


print "\nTest 7: Compute populations and densty matrices: now with phase"
el = Electronic(3,1, 0.66)
print "nstates = ", el.nstates
print "istate = ", el.istate
for i in range(0,el.nstates):
    print i, el.c(i)
    print "q[%3i]= %8.5f , p[%3i]= %8.5f, |c[%i]|^2= %8.5f" % (i, el.q[i], i, el.p[i], i, el.rho(i,i).real)

print "Density matrix:"
for i in range(0,el.nstates):
    tmp = ""
    for j in range(0,el.nstates):
        tmp = tmp + str(" %8.5f " % el.rho(i,j).real)
    print tmp






