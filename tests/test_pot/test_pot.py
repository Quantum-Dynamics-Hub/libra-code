import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
#sys.path.insert(1,cwd+"/../../_build/src/cell")
sys.path.insert(1,cwd+"/../../_build/src/pot")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygpot import *


print "\nTest 2: Importing the library and its content"
r1 = VECTOR(0.0, 0.0, 0.0)
r2 = VECTOR(1.0, 0.0, 0.0)

for i in range(0,150):
    x = 0.1 * i
    r2.x = x   
    SW, dSW = 0.0, VECTOR(0.0, 0.0, 0.0)
    SW, dSW = SWITCH(r1, r2, 5.0, 10.0 )
    print x, SW, dSW.x

