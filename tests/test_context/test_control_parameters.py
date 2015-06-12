import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/context")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygcontext import *


ctrl = ctx_Control_Parameters()
ctrl.save_xml("ctrl0.xml")


param1 = 1.0
ctrl.add(".x.param1", param1)

param2 = "Chalk"
ctrl.add(".x.param2", param2)

param3 = intList()
for i in xrange(3):
    param3.append(i)

ctrl.add(".x.param3", param3)


ctrl.save_xml("ctrl1.xml")

