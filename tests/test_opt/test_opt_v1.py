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
import cmath
import math
import os
import sys
import unittest


cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/opt")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cygopt import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libopt import *
    from liblinalg import *




def model1(q, params):
    """
    q = (x,y,z)
    E = 1/2* ((x-1)^2 + (y+2)^2 + z^2)

    """

    grd = MATRIX(3,1) 
    x,y,z = q.get(0), q.get(1), q.get(2)

    grd.set(0, x-1.0)
    grd.set(1, y+2.0)
    grd.set(2, z)

    return grd


def test():
    q = MATRIX(3,1)
    q.set(0, 1)
    q.set(1, 1)
    q.set(2, 1)

    params = {}    
           
    qopt = grad_descent(model1, q, params, 1e-6, 1.0, 10)

    print "The optimized positions are:"; qopt.show_matrix()
    

        

test()






