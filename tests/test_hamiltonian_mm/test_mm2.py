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
# Tutorial: Create ForceField object and set its parameters - also show that
# depending on functional setup you will get different number of 
# interactions
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
#sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/hamiltonian_atomistic")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygchemobjects import *
from cyghamiltonian import *
#from cyghamiltonian_atomistic import *


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
#Load_Molecule(U, syst, os.getcwd()+"/Clusters/2benz.ent", "pdb")
Load_Molecule(U, syst, os.getcwd()+"/Molecules/test1a.pdb", "pdb_1")
syst.determine_functional_groups(0)  # 

#sys1.show_atoms()
syst.show_fragments()
syst.show_molecules()

print "Number of atoms in the system = ", syst.Number_of_atoms


#======= Parameters ==============
# Create force field objects
uff = ForceField()

# Load parameters
Load_UFF(uff)

# Set up functional forms
uff.set_functionals({"bond":"Harmonic","angle":"Harmonic","vdw":"LJ12_6"})




print "\nTest 2: Object of Hamiltonian_Atomistic class"
ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)

# Choose the Hamiltonian to be of MM type
ham.set_Hamiltonian_type("MM")

# Define regions for different levels of theory
atlst1 = range(1,syst.Number_of_atoms+1)

# Create interactions between atoms from atom list 1
verb = 0
assign_rings = 0

# Set atomistically-resolved (as opposed to coarse-grained resolutions) interactions
# among all atoms in group atlst1 with all atoms in the same list
ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)    # The simplest sintax

# Show setup interactions
ham.show_interactions_statistics()




print "\nTest 3: Note that now we add dihedral functional, so eventually we will have different interactions"

uff = ForceField({"bond_functional":"Harmonic","angle_functional":"Fourier","vdw_functional":"LJ","dihedral_functional":"General0","R_vdw_on":6.0,"R_vdw_off":7.0})
Load_UFF(uff)

# Reset Hamiltonian to comply with new Force Field
ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)  # Constructing new object, otherwise the interactions will be
                                                        # added to the older one and we would overcount them
ham.set_Hamiltonian_type("MM")
ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)    # The simplest sintax
ham.show_interactions_statistics()


