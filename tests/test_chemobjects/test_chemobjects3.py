###################################################################
# Tutorial: Loading structure from files
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


# Now to Systems
i = 1
while i<=19:
    print "=====================System", i,"====================="
    # System creation section
    syst = System()

    # Create atoms, link them, define groups
    inp_file = os.getcwd()+"/Molecules/test"+str(i)+"a.pdb"
    Load_Molecule(U, syst, inp_file, "pdb_1")

    syst.determine_functional_groups(1)  # 1 - determine rings
    syst.show_atoms()
    print "======================================================="

    i = i + 1



