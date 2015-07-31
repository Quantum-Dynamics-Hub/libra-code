###################################################################
# Tutorial: Create ForceField object and set its parameters- for many systems
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


# Create force field
uff = ForceField({"bond_functional":"Harmonic","angle_functional":"Fourier",
                   "vdw_functional":"LJ","dihedral_functional":"General0",
                    "R_vdw_on":6.0,"R_vdw_off":7.0}
                )
Load_UFF(uff)
verb = 0
assign_rings = 1



#======= System ==============
for i in range(1,13):
    print "=================== System ",i,"======================="

    syst = System()
    Load_Molecule(U, syst, os.getcwd()+"/Molecules/test"+str(i)+"a.pdb", "pdb_1")

    syst.show_fragments()
    syst.show_molecules()

    print "Number of atoms in the system = ", syst.Number_of_atoms
    atlst1 = range(1,syst.Number_of_atoms+1)


    # Creating Hamiltonian
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")    
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)  
    ham.show_interactions_statistics()


