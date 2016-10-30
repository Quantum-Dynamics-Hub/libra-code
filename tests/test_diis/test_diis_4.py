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

########################################################################################
#
# In this example, we extend our line search algorithm, to minimize debug info printing
# and showcase actual working of the algorithm: finding the root of a function 
# Now we apply it to a 2D problem
# We also combine the DIIS algorithm with the gradient optimization
# To fully see the interoperation of the gradient optimization and DIIS we use more complex function
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

def func(R):
    x, y = R[0], R[1]
    f = (x - y + 0.1)**2 + 0.25*((x-1.0)**2 + y**2)**2
    dfdx =  2.0*(x-y+0.1) + ((x-1.0)**2 + y**2)*(x-1.0)
    dfdy = -2.0*(x-y+0.1) + ((x-1.0)**2 + y**2)*y

    return f, dfdx, dfdy

# In our case, f^2 is essentialy the error function

def printout(diis):
    print "DIIS coefficients", diis.get_diis_c()
    print "DIIS objective matrices and errors"; 
    X = diis.get_diis_X()
    E = diis.get_diis_err()
    sz = len(X)
    for i in xrange(sz):
        print i, X[i].get(0,0), X[i].get(1,1), E[i].get(0,0), E[i].get(1,1)


def add_set(diis, R):
    x, y = R[0], R[1]
    #print "Adding one set"
    X = MATRIX(2,2);   X.set(0,0, x);  X.set(1,1, y);  
    E = MATRIX(2,2);   z = func(R);  E.set(0,0, z[0]); E.set(1,1, z[1]); 
    diis.add_diis_matrices(X, E);
    #printout(diis)
    return z # function and its gradient


def line_search(tol, dt, R):
# This function searches for the point that delivers minimum along given direction
# The direction is given by the gradient along negative gradient of the function at
# the initial point
# The function modifies the initial R, to end up with the point that delivers minimum
# The function also returns the target function and its gradient at the final point

    diis = DIIS(3,2)
    X = MATRIX(2,2);  # extrapolated matrix
    printout(diis)


    z0 = add_set(diis, R)     # starting point 
    val0 = z0[0]

    # Initial step
    z = add_set(diis, [R[0]-dt*z0[1], R[1]-dt*z0[2]] );  val = z[0]

    diis.extrapolate_matrix(X);     
    R = [X.get(0,0), X.get(1,1)]

    err = math.fabs(val - val0);  val0 = val
    it = 0
    print it, err

    while err>tol:
        z = add_set(diis, R);  val = z[0]
        diis.extrapolate_matrix(X);    
        R = [X.get(0,0), X.get(1,1)]

        err = math.fabs(val - val0);  val0 = val

        it = it + 1
        print it, err, R

    return R, z




# This is 3 element predictor for the matrices 1x1
#diis = DIIS(3,2)
#X = MATRIX(2,2);  # extrapolated matrix
#printout(diis)

dt_ls = 0.01
dt_gr = 0.1
R = [0.0, 0.0]
tol_ls = 1e-8
tol_gr = 1e-8


it = 0
err = 2.0*tol_gr
while err>tol_gr and it<100:
    R, z = line_search(tol_ls, dt_ls, R)

    #z = func(R)
    print "Iter= ",it, " R= ", R, " z= ", z
    it = it + 1

    # Gradient step    
    sz = len(R)
    err = 0.0
    for i in xrange(sz):
        R[i] = R[i] - dt_gr * z[1+i]
        err = err + z[1+i]**2
    err = math.sqrt(err)




    

