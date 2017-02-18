#*********************************************************************************                     
#* Copyright (C) 2017 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file fit.py 
# This module implements a linear regression method and several pre-defined fit functions
#

import math

def Regression(X,Y):
# Finds y = a + b*x
    x,y,x2,y2,xy = 0.0, 0.0, 0.0, 0.0, 0.0

    sz = len(X)

    i = 0
    while i<sz:
        x = x + X[i]
        y = y + Y[i]
        x2 = x2 + X[i]*X[i]
        y2 = y2 + Y[i]*Y[i]
        xy = xy + X[i]*Y[i]
        i = i + 1
    N = float(sz)
    b = (N*xy - x*y)/(N*x2 - x*x)
    a = (y - b*x)/N

    b = y/x
    a = 0.0

    return [a,b]


def fit_exp(X,Y, x0, verbose=1):
    # Fit Y vs. X data to  Y = A*exp(-B*(X-x0))

    sz = len(X)  # The number of data points
    
    # Linearize input data
    linx = [0.0] * sz
    liny = [0.0] * sz

    for i in xrange(sz):
        linx[i] = X[i] - x0
        liny[i] = math.log(Y[i])

    # Run regression and compute parameters   
    a,b = Regression(linx,liny)  # ln(y) = ln(A) - B*(X-x0) = a + b*x, so: a = ln(A), b = -B, x = X-x0

    A = math.exp(a)
    B = -b

    # Compute prediction
    predy = [0.0] * sz
    for i in xrange(sz):
        predy[i] = A * math.exp(-B*X[i])

    if verbose:
        print "Fitting data to the expression Y = A * exp(-B*(X-x0))"
        print "With parameters: x0 = ", x0
        print "Fitting parameters: A = ", A, " B = ", B, " 1/B = ", 1.0/B

    return predy, A, B



def fit_gau(X,Y, x0, verbose=1):
    # Fit Y vs. X data to  Y = A*exp(-B*(X-x0)^2)

    sz = len(X)  # The number of data points
    
    # Linearize input data
    linx = [0.0] * sz
    liny = [0.0] * sz

    for i in xrange(sz):
        linx[i] = (X[i] - x0)**2
        liny[i] = math.log(Y[i])

    # Run regression and compute parameters   
    a,b = Regression(linx,liny)  # ln(y) = ln(A) - B*(X-x0)^2 = a + b*x, so: a = ln(A), b = -B,  x = (X-x0)^2

    A = math.exp(a)
    B = -b

    # Compute prediction
    predy = [0.0] * sz
    for i in xrange(sz):
        predy[i] = A * math.exp(-B*(X[i] - x0)**2)

    if verbose:
        print "Fitting data to the expression Y = A * exp(-B*(X-x0)**2)"
        print "With parameters: x0 = ", x0
        print "Fitting parameters: A = ", A, " B = ", B, " sqrt(B) = ", math.sqrt(B), " 1/B = ", 1.0/B, " sqrt(1.0/B) = ", math.sqrt(1.0/B)

    return predy, A, B



def get_data_from_file(filename, xindx, yindx, xminval=None, xmaxval=None, yminval=None, ymaxval=None):

    f = open(filename,"r")
    A = f.readlines()
    f.close()

    X, Y = [], []

    for a in A:
        tmp = a.split()

        x = float(tmp[xindx])
        y = float(tmp[yindx])

        is_add = 1
        if xminval != None:
            if  x < xminval:
                is_add = 0
        if xmaxval != None:
            if x > xmaxval:
                is_add = 0
        if yminval != None:
            if  y < yminval:
                is_add = 0
        if ymaxval != None:
            if y > ymaxval:
                is_add = 0

        if is_add:
            X.append(x)  
            Y.append(y)

    return X, Y


# Example of usage:
if __name__ == '__main__': 
    X, Y = get_data_from_file("relax.txt", 0, 1)
    #Ypred, A, B = fit_exp(X,Y, 0.0)  # Y = A * exp(-B*X)
    Ypred, A, B = fit_gau(X,Y, 0.0)  # Y = A * exp(-B*X^2)


    f = open("relax_and_fit.txt","w")
    for i in xrange(len(X)):
        f.write("%8.5f  %8.5f  %8.5f\n" % ( X[i],Y[i],Ypred[i]) )
    f.close()

    