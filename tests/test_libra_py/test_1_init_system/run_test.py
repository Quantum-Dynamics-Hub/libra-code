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
# "init_system" module of the "libra_py"
###################################################################

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

label = ["H"]
R = [VECTOR()]
g = [VECTOR()]
rnd = Random()
T = 300.0
sigma = 0.1
df = 1

x = init_system.init_system(label, R, g, rnd, T, sigma, df)

