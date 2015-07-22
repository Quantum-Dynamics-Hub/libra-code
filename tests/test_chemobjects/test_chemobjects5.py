###################################################################
# Tutorial: Loading files - now systems of molecules, also do grouping
###################################################################

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
from LoadPT import * # Load_PT
from LoadMolecule import * # Load_Molecule



##############################################################

# Create Universe and populate it
U = Universe()
verbose = 0
Load_PT(U, "elements.dat", verbose)


#======= Systems ==============
# System 1
print "============ Deal with system 1 (2 water molecules) =============="
syst = System()
Load_Molecule(U, syst, os.getcwd()+"/Clusters/2waters.ent", "pdb")
syst.determine_functional_groups(1)  # 1 - determine rings

#syst.show_atoms()
syst.show_fragments()
syst.show_molecules()

# System 2
print "============ Deal with system 2 (2 benzene molecules) =============="
syst1 = System()
Load_Molecule(U, syst1, os.getcwd()+"/Clusters/2benz.ent", "pdb")
syst1.determine_functional_groups(0)  # 

#syst1.show_atoms()
syst1.show_fragments()
syst1.show_molecules()





