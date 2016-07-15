#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

########################################################################################
#
# In this example, we show basic workflow of the operations with the DIIS procedure object.
# We use it for finding a root of a function. We print all the useful information to 
# give better understanding of the procedures and variables
#
########################################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


# Here we will demonstrate how to apply DIIS for efficient line search in 1D 
# so we use 1x1 matrices

def func(x):
    return (x-0.456)**2



def printout(diis):
    print "DIIS coefficients", diis.get_diis_c()
    print "DIIS objective matrices and errors"; 
    X = diis.get_diis_X()
    E = diis.get_diis_err()
    sz = len(X)
    for i in xrange(sz):
        print i, X[i].get(0), E[i].get(0)



# This is 3 element predictor for the matrices 1x1
diis = DIIS(3,1)
printout(diis)

print "Adding one set"
x = MATRIX(1,1); x.set(0, 0.0); f = MATRIX(1,1); f.set(0, func(0.0)**2); diis.add_diis_matrices(x, f);
printout(diis)

print "Adding the second set"
x = MATRIX(1,1); x.set(0, 1.0); f = MATRIX(1,1); f.set(0, func(1.0)**2); diis.add_diis_matrices(x, f);
printout(diis)

print "The extrapolated objective matrix\n"
x = MATRIX(1,1); diis.extrapolate_matrix(x); x.show_matrix()

rt = x.get(0)

for n in xrange(10):
    print "==== Cycle ", n, " ========="
    print "Adding next set"
    x = MATRIX(1,1); x.set(0, rt); f = MATRIX(1,1); f.set(0, func(rt)); diis.add_diis_matrices(x, f);
    printout(diis)

    print "The extrapolated objective matrix\n"
    x = MATRIX(1,1); diis.extrapolate_matrix(x); x.show_matrix()
    rt = x.get(0)



    

