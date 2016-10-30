#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

###################################################################
# Tutorial: Create ForceField object and set its parameters- for many systems
###################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



##############################################################

# Create Universe and populate it
U = Universe()
verbose = 0
LoadPT.Load_PT(U, "elements.dat", verbose)


# Create force field
uff = ForceField({"bond_functional":"Harmonic","angle_functional":"Fourier",
                   "vdw_functional":"LJ","dihedral_functional":"General0",
                    "R_vdw_on":6.0,"R_vdw_off":7.0}
                )
LoadUFF.Load_UFF(uff)
verb = 0
assign_rings = 1



#======= System ==============
for i in range(1,13):
    print "=================== System ",i,"======================="

    syst = System()
    LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/Molecules/test"+str(i)+"a.pdb", "pdb_1")

    syst.show_fragments()
    syst.show_molecules()

    print "Number of atoms in the system = ", syst.Number_of_atoms
    atlst1 = range(1,syst.Number_of_atoms+1)


    # Creating Hamiltonian
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")    
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)  
    ham.show_interactions_statistics()


