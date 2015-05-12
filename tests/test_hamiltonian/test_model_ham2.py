import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
#sys.path.insert(1,cwd+"/../../_build/src/cell")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cyghamiltonian import *


print "\nTest 2: Object of Hamiltonian_Model class"
ham = Hamiltonian_Model(1)

print "\nTest 3: set_params()"
ham.set_rep(1)  # adiabatic
#ham.set_params([0.1, 0.1, 1.6, 1.0])

print "\nTest 4: compute_diabatic()"
ham.compute_diabatic([1.0],[0.0])

print "\nTest 5: compute_adiabatic()"
ham.compute_adiabatic([1.0],[0.0])

print "\nTest 6: Adiabatic PESs"
f = open("dac_adia.txt","w")
for i in range(-100,100):
    x = 0.1 * i    
    ham.compute_adiabatic([x],[0.0])
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (x, ham.H(0,0).real, ham.H(1,1).real, ham.H(0,1).real ) )
#    H1 = res[2]
#    f1.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H1.get(0,0), H1.get(1,1), H1.get(0,1) ) )
#    H2 = res[3]
#    f2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H2.get(0,0), H2.get(1,1), H2.get(0,1) ) )

f.close()
#f1.close()
#f2.close()






