import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/chemobjects")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygchemobjects import *

uni1 = Universe()

at = Atom(uni1)
at.Atom_id = 1
at.globAtom_Index = 0

at.show_info()

at.save("atom.xml")


