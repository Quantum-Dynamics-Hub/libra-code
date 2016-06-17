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
# This script demonstrates the construction and use of the SD (Slater Determinant)
# data type 
###################################################################

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

# Our pool of orbitals has 4 MOs in each spin channel: e.g. HOMO-1, HOMO, LUMO, LUMO+1
mo_pool_alp = CMATRIX(4,4)
mo_pool_bet = CMATRIX(4,4)

for i in xrange(4):
    mo_pool_alp.set(i,i, 1.0, 0.0)
    mo_pool_bet.set(i,i, 1.0, 0.0)

# Consider 2 electrons in the active space:
# Ground state: 2 electrons are sitting on HOMO (HOMO-1 is not included into our active space)
#                                   alpha               beta
SD0 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1]), Py2Cpp_int([1]) )
SD1 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1]), Py2Cpp_int([2]) )

print "Ground state"
print "Normalization factor= ", SD0.normalization_factor()
print "Orbitals:\n"; SD0.get().show_matrix()
print "Spins", Cpp2Py(SD0.spin)

print "Excited state"
print "Normalization factor= ", SD1.normalization_factor()
print "Orbitals:\n"; SD1.get().show_matrix()
print "Spins", Cpp2Py(SD1.spin)


