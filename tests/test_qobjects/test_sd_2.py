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
# Here, we check the possibility of the singlet-triplet coupling
# Well, orbital mixing without breaking spatial symmetry of the alha/beta
# orbitals doesn't help
###################################################################

import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


# Our pool of orbitals has 4 MOs in each spin channel: e.g. HOMO-1, HOMO, LUMO, LUMO+1
mo_pool_alp = CMATRIX(4,4)
mo_pool_bet = CMATRIX(4,4)

# Create a situation when orbitals 1 and 2 are mixed (e.g. interacting), but the spatial
# components for alpha-channel are the same as the spatial components of the beta-channel

for i in [0,3]:
    mo_pool_alp.set(i,i, 1.0, 0.0)
    mo_pool_bet.set(i,i, 1.0, 0.0)

phi = 1.0
c, s = math.cos(phi), math.sin(phi)

mo_pool_alp.set(1,1, c, 0.0);    mo_pool_alp.set(1,2, s, 0.0);
mo_pool_alp.set(2,1, -s, 0.0);   mo_pool_alp.set(2,2, c, 0.0);

mo_pool_bet.set(1,1, c, 0.0);    mo_pool_bet.set(1,2, s, 0.0);
mo_pool_bet.set(2,1, -s, 0.0);   mo_pool_bet.set(2,2, c, 0.0);


# Consider 2 electrons in the active space:
# Ground state: 2 electrons are sitting on HOMO (HOMO-1 is not included into our active space)
# Schematics: U - spin-up, D - spin-down, I - inactive electrons
#
#       SD0            SD1           SD2
#
#  3  ---------     ---------      ---------
#  2  ---------     --   D --      -- U   --
#  1  -- U D --     -- U   --      -- U   --
#  0  -- I I --     -- I I --      -- I I --

#
#
#                                   alpha               beta
#
SD0 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1]), Py2Cpp_int([1]) ) # GS
SD1 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1]), Py2Cpp_int([2]) ) # H->L, singlet
SD2 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1,2]), Py2Cpp_int([]) ) # H->L, triplet


print "Ground state"
print "Normalization factor= ", SD0.normalization_factor()
print "Orbitals:\n"; SD0.get().show_matrix()
print "Spins =", Cpp2Py(SD0.spin)

print "Excited state, singlet"
print "Normalization factor= ", SD1.normalization_factor()
print "Orbitals:\n"; SD1.get().show_matrix()
print "Spins = ", Cpp2Py(SD1.spin)

print "Excited state, triplet"
print "Normalization factor= ", SD2.normalization_factor()
print "Orbitals:\n"; SD2.get().show_matrix()
print "Spins = ", Cpp2Py(SD2.spin)


# Now compute a overlap matrix of the above defined SDs
sd1 = SDList()
sd1.append(SD0)
sd1.append(SD1)
sd1.append(SD2)

ovlp = SD_overlap(sd1, sd1)  # overlap of the SDs at the same time, for the same (sub)system
ovlp.show_matrix()




