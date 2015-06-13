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

print "\nTest 2: The purpose of this tutorial is just to generate a template XML consistent with expected internal structure"
ctx = ctx_Control_Parameters()
print "path=", ctx.get_path()
ctx.save_xml("control_parameters.xml")
