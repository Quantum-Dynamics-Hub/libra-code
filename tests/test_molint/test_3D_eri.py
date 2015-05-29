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



print "\nTest 2: Kinetic and nuclear energy of normalized Gaussians"
f = open("3D_eri.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(1.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

for i in range(0,50):
    Rb.x = 0.1 * i

    sasa_sbsb   = electron_repulsion_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  0,0,0, 1.3, Rb)
    sasa_sasb   = electron_repulsion_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Ra, 0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb)
    sasb_sasb   = electron_repulsion_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, 0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb)

    sapxa_sbpxb = electron_repulsion_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  1,0,0, 1.3, Rb)
    sapxa_sapxb = electron_repulsion_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    sapxa_sbpxa = electron_repulsion_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  1,0,0, 1.3, Ra)
    pxapxa_sbpxa = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  1,0,0, 1.3, Ra)


    pxapxa_pxbpxb   = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 1,0,0, 1.3, Rb,  1,0,0, 1.3, Rb)
    pxapxa_pxapxb   = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    pxapxb_pxapxb   = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, 1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (Rb.x, sasa_sbsb, sasa_sasb, sasb_sasb, sapxa_sbpxb, sapxa_sapxb, sapxa_sbpxa, pxapxa_sbpxa, pxapxa_pxbpxb, pxapxa_pxapxb, pxapxb_pxapxb) )


f.close()



