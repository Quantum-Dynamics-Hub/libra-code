#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
 This file demonstrates some of the conversions between Python and C++ data types

 Therefore, one doesn't always need to pass the Pythonic objects to the Libra functions
 exposed to Python. Or, in other words, one doesn't have to expose the Libra's functions
 with the Pythonic arguments. Instead, one can define the functions with regular C++
 date types and expose them to Python.

"""


import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


x = doubleList()
x.append(1.0)
x.append(-2.0)
print x, x[0], x[1]

X = Cpp2Py(x)
print X


y = doubleList()
y.append(-1.0)
y.append( 2.0)
print y, y[0], y[1]


y1 = Py2Cpp_double([1.0, 2.0, 3.0, 4.0])
print y1, y1[0], y1[1], y1[2]


z = doubleMap()
z.append(x)
z.append(y)
print z, z[0][0], z[0][1], z[1][0], z[1][1]


a = StringDoubleMap()
a = {"A":1, "q0": 0.1}
print a["A"], a["q0"]




