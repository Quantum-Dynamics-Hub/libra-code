import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/qchem")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygqchem import *



print "\nTest 2: Pseudopotentials with 3D Gaussians (normalized)"
f = open("3D_pseudopot02.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(1.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

for i in range(0,50):
    Rb.x = 0.1 * i

    ss_pp   = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb )  # s(H)-s(H)
    spx_pp  = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # s(H)-px(H)
    spy_pp  = pseudopot02(1.0, 1.0, 2.0,Rc,   0,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # s(H)-py(H)
    pxpx_pp = pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb )  # px(H)-px(H)
    pxpy_pp = pseudopot02(1.0, 1.0, 2.0,Rc,   1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb )  # px(H)-py(H)


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (Rb.x, ss_pp, spx_pp, spy_pp, pxpx_pp, pxpy_pp ) )


f.close()



