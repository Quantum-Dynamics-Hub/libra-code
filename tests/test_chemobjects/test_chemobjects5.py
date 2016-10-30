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
# Tutorial: Loading files - now systems of molecules, also do grouping
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


#======= Systems ==============
# System 1
print "============ Deal with system 1 (2 water molecules) =============="
syst = System()
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/Clusters/2waters.ent", "pdb")
syst.determine_functional_groups(1)  # 1 - determine rings

#syst.show_atoms()
syst.show_fragments()
syst.show_molecules()

# System 2
print "============ Deal with system 2 (2 benzene molecules) =============="
syst1 = System()
LoadMolecule.Load_Molecule(U, syst1, os.getcwd()+"/Clusters/2benz.ent", "pdb")
syst1.determine_functional_groups(0)  # 

#syst1.show_atoms()
syst1.show_fragments()
syst1.show_molecules()





