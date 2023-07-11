#*********************************************************************************                     
#* Copyright (C) 2017-2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: fit
   :platform: Unix, Windows
   :synopsis: This module implements a linear regression method and several pre-defined fit functions
.. moduleauthor:: Alexey V. Akimov

"""

import math

def Regression(X,Y, opt=0):
    """

    Performs a linear regression of the Y vs. X of the form:    

        Y = b*x      opt=0  (default)
        Y = a + b*x  opt=1

    Args:
        X ( list of doubles ): The x axis data
        Y ( list of doubles ): The y axis data
        opt ( int ): Selects the type of fitting to perform:

            * opt = 0: The fitting is of the form: Y = b*X  (default)
            * opt = 1: The fitting is of a more general form: Y = a + b*X

    Returns:
        [double, double]: the linear regression coefficients a and b

            For opt==0, a = 0 

    """
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

    a,b = 0.0, 0.0
    if opt==0:
        b = y/x
        a = 0.0
    elif opt==1:
        b = (N*xy - x*y)/(N*x2 - x*x)
        a = (y - b*x)/N


    return [a,b]



def fit_exp(X,Y, x0=0.0, verbose=0, linreg_opt=0):
    """Fits Y vs. X data using:  Y = A*exp(-B*(X-x0))

    Args: 
        X ( list of doubles ): The x axis data
        Y ( list of doubles ): The y axis data
        x0 ( double ): the "origin" of the exponential function [default: 0.0]
        verbose ( int ): whether to print extra info (1) or not (0) [default: 0]
        linreg_opt ( int ): the type of the linear regression to use
            SeeAlso:: :funct:`Regression` [default: 0]

            * opt = 0: ln(Y) = b*(X-x0)
            * opt = 1: ln(Y) = a + b*(X-x0)

    Returns:
        predy ( list of doubles ): The Y values predicted using the obtained linear 
            interpolation parameters. These values are computed for all X values
   
        A ( double ):  The A coefficient in: Y = A*exp(-B*(X-x0))
        B ( double ):  The B coefficient in: Y = A*exp(-B*(X-x0))


    Example:

        >>> # get the first two columns of data as X and Y
        >>> X, Y = get_data_from_file("relax.txt", 0, 1) 
        >>>
        >>> # Y = exp(-B*X), the minimaal version
        >>> Ypred, A, B = fit_exp(X,Y)        
        >>>
        >>> # Y = exp(-B*X), the more explicit version
        >>> Ypred, A, B = fit_exp(X,Y, 0.0)   
        >>>
        >>> # Y = A * exp(-B*X), the most explicit version
        >>> Ypred, A, B = fit_exp(X,Y, 0.0, 1, 1)  
        >>>
        >>> # Assume the fitting function is: Y = A * exp(-B*(X+1.5))
        >>> Ypred, A, B = fit_exp(X,Y, 1.5, 1, 1)  # Y = A * exp(-B*(X+1.5))

        
    """

    sz = len(X)  # The number of data points
    
    # Linearize input data
    linx = [0.0] * sz
    liny = [0.0] * sz

    for i in range(0,sz):
        linx[i] = X[i] - x0
        liny[i] = math.log(Y[i])

    # Run regression and compute parameters   
    a,b = Regression(linx,liny,linreg_opt)  # ln(y) = ln(A) - B*(X-x0) = a + b*x, so: a = ln(A), b = -B, x = X-x0

    A = math.exp(a)
    B = -b

    # Compute prediction
    predy = [0.0] * sz
    for i in range(0,sz):
        predy[i] = A * math.exp(-B*X[i])

    if verbose:
        print("Fitting data to the expression Y = A * exp(-B*(X-x0))")
        print("With parameters: x0 = ", x0)
        print("Fitting parameters: A = ", A, " B = ", B, " 1/B = ", 1.0/B)

    return predy, A, B



def fit_gau(X,Y, x0=0.0, verbose=0, linreg_opt=0):
    """Fits Y vs. X data using:  Y = A*exp(-B*(X-x0)^2)

    Args: 
        X ( list of doubles ): The x axis data
        Y ( list of doubles ): The y axis data
        x0 ( double ): the "origin" of the exponential function [default: 0.0]
        verbose ( int ): whether to print extra info (1) or not (0) [default: 0]
        linreg_opt ( int ): the type of the linear regression to use
            SeeAlso:: :funct:`Regression` [default: 0]

            * opt = 0: ln(Y) = b*(X-x0)^2
            * opt = 1: ln(Y) = a + b*(X-x0)^2

    Returns:
        predy ( list of doubles ): The Y values predicted using the obtained linear 
            interpolation parameters. These values are computed for all X values
   
        A ( double ):  The A coefficient in: Y = A*exp(-B*(X-x0)^2)
        B ( double ):  The B coefficient in: Y = A*exp(-B*(X-x0)^2)

    Example:

        >>> # get the first two columns of data as X and Y
        >>> X, Y = get_data_from_file("relax.txt", 0, 1) 
        >>>
        >>> # Y = exp(-B*X^2), the minimaal version
        >>> Ypred, A, B = fit_gau(X,Y)        
        >>>
        >>> # Y = exp(-B*X^2), the more explicit version
        >>> Ypred, A, B = fit_gau(X,Y, 0.0)   
        >>>
        >>> # Y = A * exp(-B*X^2), the most explicit version
        >>> Ypred, A, B = fit_gau(X,Y, 0.0, 1, 1)  
        >>>
        >>> # Assume the fitting function is: Y = A * exp(-B*(X-1)^2)
        >>> Ypred, A, B = fit_exp(X,Y, -1.0, 1, 1)  # Y = A * exp(-B*(X-1)^2)


        
    """

    sz = len(X)  # The number of data points
    
    # Linearize input data
    linx = [0.0] * sz
    liny = [0.0] * sz

    for i in range(0,sz):
        linx[i] = (X[i] - x0)**2
        liny[i] = math.log(Y[i])

    # Run regression and compute parameters   
    a,b = Regression(linx,liny,linreg_opt)  # ln(y) = ln(A) - B*(X-x0)^2 = a + b*x, so: a = ln(A), b = -B,  x = (X-x0)^2

    A = math.exp(a)
    B = -b

    # Compute prediction
    predy = [0.0] * sz
    for i in range(0,sz):
        predy[i] = A * math.exp(-B*(X[i] - x0)**2)

    if verbose:
        print("Fitting data to the expression Y = A * exp(-B*(X-x0)**2)")
        print("With parameters: x0 = ", x0)
        print("Fitting parameters: A = ", A, " B = ", B, " sqrt(B) = ", math.sqrt(B), " 1/B = ", 1.0/B, " sqrt(1.0/B) = ", math.sqrt(1.0/B) )

    return predy, A, B


# Example of usage:
if __name__ == '__main__': 
    X, Y = get_data_from_file("relax.txt", 0, 1)
    Ypred, A, B = fit_exp(X,Y, 0.0)  # Y = exp(-B*X)
    #Ypred, A, B = fit_gau(X,Y, 0.0)  # Y = exp(-B*X^2)

    #Ypred, A, B = fit_exp(X,Y, 0.0, 1, 1)  # Y = A * exp(-B*X)
    #Ypred, A, B = fit_gau(X,Y, 0.0, 1, 1)  # Y = A * exp(-B*X^2)

    #Ypred, A, B = fit_exp(X,Y, 1.5, 1, 1)  # Y = A * exp(-B*(X+1.5))
    #Ypred, A, B = fit_gau(X,Y,-1.0, 1, 1)  # Y = A * exp(-B*(X-1.0)^2)




    f = open("relax_and_fit.txt","w")
    for i in range(0,len(X)):
        f.write("%8.5f  %8.5f  %8.5f\n" % ( X[i],Y[i],Ypred[i]) )
    f.close()

    