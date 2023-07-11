#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import cmath
import math
import os
import sys
import unittest
import random

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")
sys.path.insert(1,cwd+"/../../_build/src/basis")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cyglinalg import *
    from cygbasis import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from liblinalg import *
    from libbasis import *


basis = basisset_struct(54, 54, 100);

#sys.exit(0)

basis.read_basisset_file("STO-3G", 0)

#sys.exit(0)

for i in xrange(3):
    print "==== atom ", i, " ======"
    print basis.atoms[i]
    print "# of shells ", basis.atoms[i].noOfShells
    print "# of shells ", basis.atoms[i].shells
    print len(basis.atoms[i].shells)

    for j in xrange(1):
        print "shell ", basis.atoms[i].shells[j]
        print "type = ", basis.atoms[i].shells[j].shell_type
        print "contrCount = ", basis.atoms[i].shells[j].contrCount
        print "shell_ID = ", basis.atoms[i].shells[j].shell_ID

        for k in xrange(basis.atoms[i].shells[j].contrCount):
            print "Exp = ", basis.atoms[i].shells[j].exponentList[k]
            print "Coeff = ", basis.atoms[i].shells[j].coeffList[k]




