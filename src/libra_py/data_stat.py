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

import common_utils as comn




def vec_average(data):
    """This function computes the average value of the data series

    Args:
        data ( list of VECTOR objects ): sequence of real-valued ndof-dimensional vectors

    Returns:
        VECTOR: the average value for each DOF in the time-series
    """

    ave = VECTOR()
    sz = len(data)

    for i in xrange(sz):
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
    for i in xrange(sz):
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
    for i in xrange(1, sz-1):
        if data[i] > data[i-1] and data[i] > data[i+1]:
            max_indxs.append(i)

    inp = []
    sz = len(max_indxs)
    for i in xrange(sz):
        inp.append( [ max_indxs[i], data[max_indxs[i]] ] )

    out = merge_sort(inp)  # largest in the end

    if do_output:
        lgfile = open(logname, "a")
        lgfile.write("Found maxima of the data:\n")
        for i in xrange(sz):
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
    for i in xrange(N):
        res = res + data[i]
    res = res / float(N)

    #===== Std ========
    res2 = 0.0
    
    for i in xrange(N):
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
    for i in xrange(N):
        res = res + X[i]
    res = res / float(N)


    #===== Std ========
    res2 = MATRIX(X[0]); res2 *= 0.0

    for a in xrange(res2.num_of_rows):
        for b in xrange(res2.num_of_cols):
        
            tmp = 0.0
            for i in xrange(N):
                tmp = tmp + (X[i].get(a,b) - res.get(a,b))**2
            tmp = math.sqrt( tmp / float(N) )

            res2.set(a,b, tmp)

    # Find maximal and minimal values
    up_bound = MATRIX(X[0]); up_bound *= 0.0
    dw_bound = MATRIX(X[0]); dw_bound *= 0.0


    for a in xrange(res2.num_of_rows):
        for b in xrange(res2.num_of_cols):

            up_bound.set(a,b, X[0].get(a,b))
            dw_bound.set(a,b, X[0].get(a,b))

            for i in xrange(N):
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
    for i in xrange(N):
        res = res + X[i]
    res = res / float(N)


    #===== Std ========
    res2 = CMATRIX(X[0]); res2 *= 0.0

    for a in xrange(res2.num_of_rows):
        for b in xrange(res2.num_of_cols):
        
            tmp = 0.0+0.0j
            for i in xrange(N):
                dx = X[i].get(a,b) - res.get(a,b)
                tmp = tmp + (dx.conjugate() * dx)
            tmp = math.sqrt( tmp / float(N) )

            res2.set(a,b, tmp)

    return res, res2



