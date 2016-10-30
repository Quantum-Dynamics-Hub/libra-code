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
# Tutorial: Create ForceField object and set its parameters
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


#======= System ==============
syst = System()
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/Clusters/2benz.ent", "pdb")
syst.determine_functional_groups(0)  # 

#sys1.show_atoms()
syst.show_fragments()
syst.show_molecules()



#======= Parameters ==============
# Create force field objects
uff = ForceField()
# Load parameters
LoadUFF.Load_UFF(uff)

# Set up functional forms
uff.set_functionals({"bond":"Harmonic","angle":"Harmonic","vdw":"LJ12_6"})


