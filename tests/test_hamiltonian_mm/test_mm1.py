#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

###################################################################
# Tutorial: Create ForceField object and set its parameters
###################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/chemobjects")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")


print "\nTest 1: Importing the library and its content"
from cygmmath import *

from cygchemobjects import *
from LoadPT import * # Load_PT
from LoadMolecule import * # Load_Molecule
from LoadUFF import*



##############################################################

# Create Universe and populate it
U = Universe()
verbose = 0
Load_PT(U, "elements.dat", verbose)


#======= System ==============
syst = System()
Load_Molecule(U, syst, os.getcwd()+"/Clusters/2benz.ent", "pdb")
syst.determine_functional_groups(0)  # 

#sys1.show_atoms()
syst.show_fragments()
syst.show_molecules()



#======= Parameters ==============
# Create force field objects
uff = ForceField()
# Load parameters
Load_UFF(uff)

# Set up functional forms
uff.set_functionals({"bond":"Harmonic","angle":"Harmonic","vdw":"LJ12_6"})


