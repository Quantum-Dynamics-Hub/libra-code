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

##
# \file test.py This file demonstrates how to read QE wavefunctions using
# Context class of the Libra package. We also show how to use the created objects
#
#  --- No SOC case --- 
# Here, we also will consider only those orbitals that we really need
# 
# This example can also be found in test_qe directory - see the examples
# only in those cases, all the internals of the functions are determined in the 
# example script explicitly.
# Here, we test those  functions adopted in the main Libra package
#

import os
import math

# Fisrt, we add the location of the library to test to the PYTHON path
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



print "Reading energies"
e = QE_methods.read_qe_index("x.export/index.xml", [6,7] )
e.show_matrix()


coeff_1 = QE_methods.read_qe_wfc("x.export/wfc.1", "Kpoint.1", [6,7] )
coeff_2 = QE_methods.read_qe_wfc("x.export/wfc.2", "Kpoint.2", [6,7] )

nbnd = coeff_1.num_of_cols
print "The number of bands = ", nbnd
 
ovlp_1  = coeff_1.H() * coeff_1
ovlp_2  = coeff_2.H() * coeff_2
ovlp_12 = coeff_1.H() * coeff_2

print "<alp|alp> = "; ovlp_1.show_matrix()
print "<bet|bet> = "; ovlp_2.show_matrix()
print "<alp|bet> = "; ovlp_12.show_matrix()
print "\n"


for n in xrange(nbnd):
    print n, ovlp_1.get(n,n), ovlp_2.get(n,n), ovlp_12.get(n,n)  #  ovlp.get(2*n,2*n) + ovlp.get(2*n+1,2*n+1)





