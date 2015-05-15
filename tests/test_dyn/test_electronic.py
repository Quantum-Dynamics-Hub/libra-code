import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/dyn")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygdyn import *


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


print "\nTest 5: Compute populations and densty matrices"
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







