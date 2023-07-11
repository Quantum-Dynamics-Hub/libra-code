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

######################################################################
#
# Here we try to work on some acceleration of computations
#
######################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



for x in [40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 200.0, 500.0, 1000.0 ]:
    print x, Fn(0, x), math.exp(-x)

#print Fn(0, 70.0)


