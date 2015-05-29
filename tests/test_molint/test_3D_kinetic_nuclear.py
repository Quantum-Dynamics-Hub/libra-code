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
f = open("3D_kin_nucl.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(1.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

for i in range(0,50):
    Rb.x = 0.1 * i

    s_s   = kinetic_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb)
    s_s   = s_s + nuclear_attraction_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, Ra)
    s_s   = s_s + nuclear_attraction_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, Rb)

    s_px   = kinetic_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    s_px   = s_px + nuclear_attraction_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Ra)
    s_px   = s_px + nuclear_attraction_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Rb)

    px_px   = kinetic_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    px_px   = px_px + nuclear_attraction_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Ra)
    px_px   = px_px + nuclear_attraction_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Rb)

    px_py   = kinetic_integral(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb)
    px_py   = px_py + nuclear_attraction_integral(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, Ra)
    px_py   = px_py + nuclear_attraction_integral(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, Rb)

    f_f   = kinetic_integral(3,0,0, 1.3, Ra,  3,0,0, 1.3, Rb)
    f_f   = f_f + nuclear_attraction_integral(3,0,0, 1.3, Ra,  3,0,0, 1.3, Rb, Ra)
    f_f   = f_f + nuclear_attraction_integral(3,0,0, 1.3, Ra,  3,0,0, 1.3, Rb, Rb)


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (Rb.x, s_s, s_px, px_px, px_py, f_f ) )


f.close()



