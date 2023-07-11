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
# Tutorial: This file demonstrates the functionality of the 
# "init_ensembles" module of the "libra_py"
###################################################################

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


print "============= Step 1 =============="
label = ["H"]
R = [VECTOR()]
g = [VECTOR()]
rnd = Random()
T = 300.0
sigma = 0.1
df = 1

x1 = init_system.init_system(label, R, g, rnd, T, sigma, df)
x2 = init_system.init_system(label, R, g, rnd, T, sigma, df)

print "============= Step 2 ==============="

ntraj = 2
nnucl = 3*1
nel = 2
verbose = 1


print "One version"
ham, ham_adi, d1ham_adi, ham_vib = init_ensembles.init_ext_hamiltonians(ntraj, nnucl, nel)

print "The version with printing"
ham, ham_adi, d1ham_adi, ham_vib = init_ensembles.init_ext_hamiltonians(ntraj, nnucl, nel, verbose)


print "============= Step 3 ==============="

mol = init_ensembles.init_mols([x1, x2], ntraj, nnucl, verbose)


print "============= Step 4 ==============="

params = { "nu_therm": 0.01,
           "NHC_size": 3, 
           "Temperature": 300.0, 
           "thermostat_type": "Nose-Hoover"
         }

therms = init_ensembles.init_therms(ntraj, nnucl, params, verbose)


