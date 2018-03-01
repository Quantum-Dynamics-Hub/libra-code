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
sys.path.insert(1,cwd+"/../../_build/src/montecarlo")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")
sys.path.insert(1,cwd+"/../../_build/src/math_random")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cygmontecarlo import *
    from cyglinalg import *
    from cygrandom import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libmontecarlo import *
    from liblinalg import *
    from librandom import *




def piab(q, params):
    """
    The probability density function

    """

    L = params["L"]
    n = params["n"]

    x = q.get(0)

    p = 0.0
    if x>0.0 and x<L:
        p = (math.sin(x*n*math.pi/L))**2


    return p


def bin(sample, min_, max_, dx, i, j, filename):

    # Prepare the grids
    x_points, y_points = [], []
    max_pts = int((max_ - min_)/dx) + 1

    for n in xrange(max_pts):
        x_points.append(min_ + n * dx)
        y_points.append(0.0)

    # Compute the frequencies
    sz = len(sample)
    for n in xrange(sz):
        x = sample[n].get(i,j) 
        indx = int((x - min_)/dx)

        y_points[indx] = y_points[indx] + 1.0/float(sz)
              

    f = open(filename, "w")
    for n in xrange(max_pts):
        f.write("%8.5f  %8.5f \n" % (x_points[n], y_points[n]))
    f.close()


def test():

    rnd = Random()
    params = {"L":1.0, "n":5}    
    q = MATRIX(1,1);  q.set(0, 0.5)
           
    sample = metropolis_gau(rnd, piab, q, params, 1000000, 10, 0.05) 
    print len(sample)
    bin(sample, -1.0, 2.0, 0.01, 0, 0, "_distrib.txt")    


           

test()






