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
# In this example, we extend our line search algorithm, to minimize debug info printing
# and showcase actual working of the algorithm: finding the root of a function 
# Now we apply it to a 2D problem - it is essential to run the algorithm along the gradient direction
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

def func(x,y):
    return (x-0.456)**2 + (y+0.56)**2, 2.0*(x-0.456), 2.0*(y+0.56)

# In our case, f^2 is essentialy the error function

def printout(diis):
    print "DIIS coefficients", diis.get_diis_c()
    print "DIIS objective matrices and errors"; 
    X = diis.get_diis_X()
    E = diis.get_diis_err()
    sz = len(X)
    for i in xrange(sz):
        print i, X[i].get(0,0), X[i].get(1,1), E[i].get(0,0), E[i].get(1,1)


def add_set(x,y):
    #print "Adding one set"
    X = MATRIX(2,2);   X.set(0,0, x);  X.set(1,1, y);  
    E = MATRIX(2,2);   z = func(x,y);  E.set(0,0, z[0]); E.set(1,1, z[1]); 
    diis.add_diis_matrices(X, E);
    #printout(diis)
    return z # function and its gradient



# This is 3 element predictor for the matrices 1x1
diis = DIIS(3,2)
printout(diis)

dt = 1.0
#z0 = add_set( -1.0, -1.0)  # too good starting point
z0 = add_set( 0.0, 0.0)     # less efficient starting point

# add few steps along the direction of the negative gradient
add_set( -1.0-dt*z0[1], -1.0-dt*z0[2])
printout(diis)

add_set( -1.0-2.0*dt*z0[1], -1.0-2.0*dt*z0[2])
printout(diis)


print "The extrapolated objective matrix\n"
X = MATRIX(2,2); diis.extrapolate_matrix(X); X.show_matrix()

rt = [X.get(0,0), X.get(1,1)]

# Line search along the negative gradient direction
for n in xrange(25):
 
    add_set(rt[0],rt[1])   
    diis.extrapolate_matrix(X); 
   
    rt = [X.get(0,0), X.get(1,1)]

    print n,rt, func(rt[0],rt[1])[0]



    

