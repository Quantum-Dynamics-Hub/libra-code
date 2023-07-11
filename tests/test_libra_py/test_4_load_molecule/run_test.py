#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
###################################################################
# Tutorial: This file demonstrates how to load molecular structures using
# the "LoadMolecule" module of the "libra_py"
###################################################################

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


# Create Universe and populate it
U = Universe()
verbose = 0
LoadPT.Load_PT(U, "elements.dat", verbose)


#======= System ==============
syst = System()

LoadMolecule.Load_Molecule(U, syst, "Clusters/2benz.ent", "pdb")

print "Show properties"
syst.show_atoms()
syst.show_fragments()
syst.show_molecules()

