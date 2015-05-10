import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
#sys.path.insert(1,cwd+"/../../_build/src/cell")
sys.path.insert(1,cwd+"/../../_build/src/operators")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygoperators import *


print "\nTest 2: Rotate"
r = [1.0, 0.0]
dphi = math.radians(10.0)
for i in range(0,10):
    r = rotate(r[0],r[1], dphi)
    phi = math.degrees(math.atan2(r[1],r[0]))
    print r, phi

print "\nTest 3: Shift"
x = 0.0
for i in range(0, 10):
    print i, x
    x = shift(x, 0.1)


print "\nTest 4: Scale"
x = 1.0
for i in range(0, 10):
    print i, x
    x = scale(x, -0.99)

    
#void shift(double&, double);
#double expt_shift(double x,double phi);

#void scale(double&, double);
#double expt_scale(double x,double phi);

