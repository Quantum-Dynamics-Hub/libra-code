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

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


print "\nTest 2: Constructing"
print "x = DATA([1.0, 0.5, 2.0, -0.5])"
x = DATA([1.0, 0.5, 2.0, -0.5])

print "\nTest 3: Estimators"
ave, var, sd, se, mse, mae, rmse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
#print "x.Calculate_Estimators(ave, var, sd, se, mse, mae, rmse)"
#x.Calculate_Estimators(ave, var, sd, se, mse, mae, rmse)
print "ave = ",ave
print "var = ",var
print "sd = ",sd
print "se = ",se
print "mse = ",mse
print "mae = ",mae
print "rmse = ",rmse




