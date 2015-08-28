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
# Tutorial: Loading files and ring determination
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


t = Timer()

# Now to Systems
i = 1
while i<=12:
    print "=====================System", i,"====================="
    # System creation section
    syst = System()

    # Create atoms, link them, define groups
    inp_file = os.getcwd()+"/Rings/test"+str(i)+".pdb"
    Load_Molecule(U, syst, inp_file, "pdb_1")



    t.start()
    syst.determine_functional_groups(1)  # 1 - determine rings
    t.stop()
    print "Time to compute = ", t.show(), " sec"

    syst.show_rings()
    syst.show_atoms()
    print "======================================================="

    i = i + 1

# 12 (Fullerene) takes 2-3 min. to determine rings, but it is doable


