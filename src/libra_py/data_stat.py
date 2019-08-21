#*********************************************************************************                     
#* Copyright (C) 2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: data_stat
   :platform: Unix, Windows
   :synopsis: 
       This module implements various functions for data analysis

.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import common_utils as comn
import util.libutil as comn




def vec_average(data):
    """This function computes the average value of the data series

    Args:
        data ( list of VECTOR objects ): sequence of real-valued ndof-dimensional vectors

    Returns:
        VECTOR: the average value for each DOF in the time-series
    """

    ave = VECTOR()
    sz = len(data)

    for i in range(0,sz):
        ave = ave + data[i] 

    ave = ave/float(sz)
    return ave


def vec_center_data(data):
    """

    This function centers data on zero, by subtracting the average 
    value from each element, dof by dof

    Args:
        data ( list of VECTOR objects ): sequence of real-valued ndof-dimensional vectors

    Returns:
        VECTOR: the fluctuation (deviation) value for each DOF in the time-series (dX(t) = X(t) -<X> )
    """

    data_new = []    
    ave = vec_average(data)
    for d in data:
        data_new.append(d - ave)

    return data_new



def mat_average(data):
    """This function computes the average value of the data series

    Args:
        data ( list of MATRIX(ndof, 1) objects ): sequence of real-valued ndof-dimensional vectors

    Returns:
        MATRIX(ndof, 1): the average value for each DOF in the time-series
    """

    ndof = data[0].num_of_rows  
    ave = MATRIX(ndof, 1)

    sz = len(data)
    for i in range(0,sz):
        ave = ave + data[i] 

    ave = ave/float(sz)
    return ave


def mat_center_data(data):
    """

    This function centers data on zero, by subtracting the average 
    value from each element, dof by dof

    Args:
        data ( list of MATRIX(ndof, 1) objects ): sequence of real-valued ndof-dimensional vectors

    Returns:
        MATRIX(ndof, 1): the fluctuation (deviation) value for each DOF in the time-series (dX(t) = X(t) -<X> )
    """

    data_new = []    
    ave = mat_average(data)
    for d in data:
        data_new.append(d - ave)

    return data_new




def find_maxima(data, params):
    """

    This function finds all the maxima of the data set and sorts them according to the data
    The maxima are defined as data[i-1] < data[i] > data[i+1]
    
    Args: 
        data ( list of doubles ): data to be analyzed
        params ( Python dictionary ): parameters controlling the execution
 
            * **params["do_output"] ( Boolean ): wheather to output the data to file [ default: False ]
            * **params["logname"] ( Boolean ): the name of the output file [ default: "data_maxima.txt" ]
            * **params["verbose"] ( int ): the amount of the debug/descriptive info to print out [ default: 0 ]

    Returns:
        list of 2-element lists: our[i][v], where:
 
            out[i][0] - containing indices of the maximal values

    Examples:

        >>> res = find_maxima( [0.01, 1.0, 0.1, 0.25, 0.5, 0.75, -1.0 ], {} )
        >>> print res
        >>> [ [1, 1.0], [5, 0.75] ]

        >>> res = find_maxima( [0.01, 0.75, 0.1, 0.25, 0.5, 1.75, -1.0 ], {} )
        >>> print res
        >>> [ [5, 1.75], [1, 0.75] ]

    """

    critical_params = []
    default_params = { "do_output":False, "logname":"data_maxima.txt", "verbose":0 }
    comn.check_input(params, default_params, critical_params)

    do_output = params["do_output"]
    logname = params["logname"]
    verbose = params["verbose"]


    max_indxs = []
    sz = len(data)
    for i in range(1, sz-1):
        if data[i] > data[i-1] and data[i] > data[i+1]:
            max_indxs.append(i)

    inp = []
    sz = len(max_indxs)
    for i in range(0,sz):
        inp.append( [ max_indxs[i], data[max_indxs[i]] ] )

    out = merge_sort(inp)  # largest in the end

    if do_output:
        lgfile = open(logname, "a")
        lgfile.write("Found maxima of the data:\n")
        for i in range(0,sz):
            lgfile.write("maximum index = %3i  index of the datapoint = %8.5f  data value = %8.5f \n" % (i, out[sz-1-i][0], out[sz-1-i][1]) )
        lgfile.close()    
    
    return out
    




def scalar_stat(data):
    """

    The function computes some simple descriptive statistics of the scalar data series.

    Args:
        data ( list of doubles ): the data to be analyzed

    Returns:
        tuple: (res, res2), where:

            * res ( double ): average of data
            * res2 ( double ): standard deviation of data

    """

    N = len(data)
    res = 0.0

    #===== Average ====
    for i in range(0,N):
        res = res + data[i]
    res = res / float(N)

    #===== Std ========
    res2 = 0.0
    
    for i in range(0,N):
        res2 = res2 + (data[i] - res)**2
    res2 = math.sqrt( res2 / float(N) )

    return res, res2




def mat_stat(X):
    """Computes the average and standard deviation of list of MATRIX(N,N) objects

    Args:
        X ( list of MATRIX(N,N) objects ): the data to be analyzed

    Returns:
        tuple: (res, res2, dw_bound, up_bound): where

            * res ( MATRIX(N,N) ): average of each matrix element of data
            * res2 ( MATRIX(N,N) ): standard deviation of each matrix element of data
            * dw_bound ( MATRIX(N,N) ): lower bounds of the data for each matrix element
            * up_bound ( MATRIX(N,N) ): upper bounds of the data for each matrix element

    """

    N = len(X)
    res = MATRIX(X[0]); res *= 0.0

    #===== Average ====
    for i in range(0,N):
        res = res + X[i]
    res = res / float(N)


    #===== Std ========
    res2 = MATRIX(X[0]); res2 *= 0.0

    for a in range(0,res2.num_of_rows):
        for b in range(0,res2.num_of_cols):
        
            tmp = 0.0
            for i in range(0,N):
                tmp = tmp + (X[i].get(a,b) - res.get(a,b))**2
            tmp = math.sqrt( tmp / float(N) )

            res2.set(a,b, tmp)

    # Find maximal and minimal values
    up_bound = MATRIX(X[0]); up_bound *= 0.0
    dw_bound = MATRIX(X[0]); dw_bound *= 0.0


    for a in range(0,res2.num_of_rows):
        for b in range(0,res2.num_of_cols):

            up_bound.set(a,b, X[0].get(a,b))
            dw_bound.set(a,b, X[0].get(a,b))

            for i in range(0,N):
                xab = X[i].get(a,b)
                if xab > up_bound.get(a,b):
                    up_bound.set(a,b, xab)
                if xab < dw_bound.get(a,b):
                    dw_bound.set(a,b,xab)


    return res, res2, dw_bound, up_bound




def cmat_stat(X):
    """Computes the average and standard deviation of list of CMATRIX(N,N) objects

    Args:
        X ( list of CMATRIX(N,N) objects ): the data to be analyzed

    Returns:
        tuple: (res, res2): where

            * res ( CMATRIX(N,N) ): average of each matrix element of data
            * res2 ( CMATRIX(N,N) ): standard deviation of each matrix element of data

    """

    N = len(X)
    res = CMATRIX(X[0]); res *= 0.0

    #===== Average ====
    for i in range(0,N):
        res = res + X[i]
    res = res / float(N)


    #===== Std ========
    res2 = CMATRIX(X[0]); res2 *= 0.0

    for a in range(0,res2.num_of_rows):
        for b in range(0,res2.num_of_cols):
        
            tmp = 0.0+0.0j
            for i in range(0,N):
                dx = X[i].get(a,b) - res.get(a,b)
                tmp = tmp + (dx.conjugate() * dx)
            tmp = math.sqrt( tmp / float(N) )

            res2.set(a,b, tmp)

    return res, res2






def cmat_stat2(X, opt):
    """Computes the norm-N average of a list of CMATRIX(N,N) objects

    Args:
        X ( list of CMATRIX(N,N) objects ): the data to be analyzed
        opt ( int ): the option for averaging:

            * opt == 0 :  t_ij = <x_ij>  + i <y_ij>
            * opt == 1 :  t_ij = <|x_ij|>  + i <|y_ij|>
            * opt == 2 :  t_ij = sqrt(<x_ij^2>)  + i sqrt(<y_ij^2>)
            * opt == 3 :  t_ij = sqrt(<x_ij^2> + <y_ij^2>) = sqrt(|z_ij|^2) 
            
    Returns:
        CMATRIX(N,N): norm-N average of each matrix element of data


    """

    ndata = len(X)
    N = X[0].num_of_cols

    res = CMATRIX(N, N)


    #===== Average ====
    for idata in range(0,ndata):

        for i in range(0,N):
            for j in range(0,N):
 
                if opt == 0:
                    res.add(i,j, X[idata].get(i,j))

                elif opt == 1:
                    re = math.fabs(X[idata].get(i,j).real)
                    im = math.fabs(X[idata].get(i,j).imag)
                    res.add(i,j, re+1.0j*im)

                elif opt == 2:
                    re = math.fabs(X[idata].get(i,j).real)
                    im = math.fabs(X[idata].get(i,j).imag)
                    res.add(i,j, re*re+1.0j*im*im)

                elif opt == 3:
                    re = math.fabs(X[idata].get(i,j).real)
                    im = math.fabs(X[idata].get(i,j).imag)
                    res.add(i,j, math.sqrt(re*re+im*im)*(1.0+0.0j) )


    res = res / float(ndata)

    for i in range(0,N):
        for j in range(0,N):

            if opt == 2:
                re = math.sqrt( res.get(i,j).real )
                im = math.sqrt( res.get(i,j).imag )
                res.set(i,j, re+1.0j*im )


    return res



def cmat_distrib(X, i, j, component, xmin, xmax, dx):
    """Computes the distribution of the matrix element values in the list of CMATRIX(N,N) objects

    Args:
        X ( list of CMATRIX(N,N) objects ): the data to be analyzed
        i ( int ): row index of the matrix element to analyze
        j ( int ): column index of the matrix element to analyze
        component ( int = 0 or 1 ): determines whether to analyze the real (0) or imaginary (1)
            component of the data series
        xmin ( double ): the minimal value of the bin support
        xmax ( double ): the maximal value of the bin support
        dx ( double ): the value of the bin support grid spacing
                    
    Returns:
        tuple: ( bin_support, dens, cum ), where:

            * bin_support ( list of doubles ): the range of values the distribution is computed for
            * dens ( list of doubles ): the probability density
            * cum ( list of doubles ): the cumulative distribution function 

    """

    #============= Extract the data ===========    
    ndata = len(X)
    y = []
    for idata in range(0,ndata):
        if component==0:
            y.append( X[idata].get(i,j).real )
        elif component==1:
            y.append( X[idata].get(i,j).imag )
    data = DATA(y)    

    #============= Build the bin support ===========
    bin_support = []
    x = xmin
    while x <= xmax:        
        bin_support.append(x)
        x = x + dx
        

    dens, cum = data.Calculate_Distribution(bin_support)

    return bin_support, dens, cum



def compute_density(X, Y, minx, maxx, dx):
    """Computes the probability density from the non-uniform distribution of pair points
   
    Args:
        X ( list of double ): the original values of x
        Y ( list of double ): the original values of y
        minx ( double ): the minimal value of the new X axis
        maxx ( double ): the max value of the new X axis
        dx ( double ): the grid spacing of the new X axis

    Note:
        The pair relationship X[i] - Y[i] is expected

    Returns:
        tuple: ( nX, nY ): where
 
            * nX ( list of doubles ) - new, uniform X axis
            * nY ( list of doubles ) - renormalized Y axis
    """

    # Prepare the grids
    nX, nY, cnt = [], [], []
    max_pts = int((maxx - minx)/dx) + 1

    for n in range(0,max_pts):
        nX.append(minx + n * dx)
        nY.append(0.0)
        cnt.append(0.0)

    # Compute the frequencies
    sz = len(X)
    for n in range(0,sz):
        indx = int((X[n] - minx)/dx)

        if indx>0 and indx<max_pts:
            nY[indx] = nY[indx] + Y[n]
            cnt[indx] = cnt[indx] + 1.0

    for n in range(0,max_pts):
        if cnt[n]>0.0:
            nY[n] = nY[n]/cnt[n]
              
    return nX, nY



def bootstrapping(X, nB, rnd):
    """Bootstrapping the sample

    Args:
        X ( list of N int/double ): the original sample

        nB ( int ): the number of bootstrapping cycles

    Returns: 
    
    """

    N = len(X) 

    # Actual re-sampling
    sample = MATRIX(N, nB)
    for samp in range(0,nB):
        for n in range(0,N):
            i = int(rnd.uniform(0, N))
            sample.set(n, samp, X[i])

    # Average over 
    sample = MATRIX(N, nB)
    





