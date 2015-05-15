import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygdyn import *
from cyghamiltonian import *


print "\nTest 2: Compute populations and densty matrices"
el = Electronic(2,0)
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


# Initialize Hamiltonian: SAC in adiabatic representation
#ham = Hamiltonian() 
x = -10.0
v = 1.0
ham = Hamiltonian_Model(0)
ham.set_rep(1)  # adiabatic
ham.set_q([x])  
ham.set_v([v])  
ham.compute()

# CPA-style dynamics: the trajectory is precomputed (although velocity is fixed)
dt = 0.1
t = 0.0
f = open("dyn_sac.txt","w")
for i in range(0,200):  # time steps
    t = t + dt
    x = x + v * dt
    ham.set_q([x])
    ham.compute()

    # Here is the electronic propagation:
    el.propagate_electronic(dt, ham)

    # Adiabatic PES and couplings
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, x, el.rho(0,0).real, el.rho(1,1).real, el.rho(0,1).real, el.rho(0,1).imag ) )

f.close()




print "\nTest 3: Rabi2 model"
el = Electronic(2,0)
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


# Initialize Hamiltonian: Rabi2 in adiabatic representation
#ham = Hamiltonian() 
x = -10.0
v = 10.0
ham = Hamiltonian_Model(5)
ham.set_params([0.1, 0.1, 1.0])
ham.set_rep(0)  # diabatic
ham.set_q([x])  
ham.set_v([v])  
ham.compute()

# CPA-style dynamics: the trajectory is precomputed (although velocity is fixed)
dt = 0.1
t = 0.0
f = open("dyn_rabi2.txt","w")
for i in range(0,200):  # time steps
    t = t + dt
    x = x + v * dt
    ham.set_q([x])
    ham.compute()

    # Here is the electronic propagation:
    el.propagate_electronic(dt, ham)

    # Adiabatic PES and couplings
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, x, el.rho(0,0).real, el.rho(1,1).real, el.rho(0,1).real, el.rho(0,1).imag ) )

f.close()

