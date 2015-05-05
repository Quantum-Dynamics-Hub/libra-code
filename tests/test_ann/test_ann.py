import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/ann")


print "\nTest 1: Importing the library and its content"
print "from cygann import *"
from cygann import *

print "\nTest 2: Constructing"
