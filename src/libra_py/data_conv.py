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
.. module:: data_conv
   :platform: Unix, Windows
   :synopsis: 
       This module implements various functions for data conversions and data transformations

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


def transform_data(X, params):
    """

    This is an auxiliary function to transform the original matrices X (e.g. H_vib) according to:     
    X(original) ->    ( X + shift1 ) (x) (scale) + shift2,  

    Here, (x) indicates the element-wise multiplicaiton, and shift1, shift2, and scale are matrices

    Args: 
        X ( list of lists of CMATRIX(nstates, nstates) ): the original data stored as
            X[idata][step] - a CMATRIX(nstates, nstates) for the dataset `idata` and time 
            step `step`
        params ( dictionary ): parameters controlling the transformation

            * **params["shift1"]** ( CMATRIX(nstates,nstates) ): first shift corrections [units of X], [default: zero]
            * **params["shift2"]** ( CMATRIX(nstates,nstates) ): second shift corrections [units of X], [default: zero]
            * **params["scale"]** ( CMATRIX(nstates,nstates) ): scaling of the X [unitless], [default: 1.0 in all matrix elements]

    Returns:    
        None: but transformes X input directly, so changes the original input


    Example:
        Lets say we have a data set of 2x2 matrices and we want to increase the energy gap by 0.1 units and scale the
        couplings by a factor of 3. Then, the input is going to be like this:

        >>> scl = CMATRIX(2,2); 
        >>> scl.set(0,0, 1.0+0.0j);  scl.set(0,1, 3.0+0.0j);
        >>> scl.set(0,0, 3.0+0.0j);  scl.set(0,1, 1.0+0.0j);
        >>> shi = CMATRIX(2,2); 
        >>> shi.set(0,0, 0.0+0.0j);  shi.set(0,1, 0.0+0.0j);
        >>> shi.set(0,0, 0.0+0.0j);  shi.set(0,1, 0.1+0.0j);
        >>> transform_data(X, {"shift2":shi, "scale":scl })


    """

    ndata = len(X)
    nsteps = len(X[0])
    nstates = X[0][0].num_of_cols

    sh1 = CMATRIX(nstates, nstates) # zero
    sh2 = CMATRIX(nstates, nstates) # zero
    scl = CMATRIX(nstates, nstates) # all elements are 1

    for i in xrange(nstates):
        for j in xrange(nstates):
            scl.set(i,j, 1.0+0.0j)

    critical_params = [  ] 
    default_params = { "shift1":sh1, "shift2":sh2, "scale":scl  }
    comn.check_input(params, default_params, critical_params)

    for idata in xrange(ndata):
        for istep in xrange(nsteps):

            tmp = CMATRIX(X[idata][istep])
            tmp = tmp + params["shift1"]
            tmp.dot_product( tmp, params["scale"] )
            tmp = tmp + params["shift2"]
            X[idata][istep] = CMATRIX(tmp)




def unpack1(H, i, j, component):
    """

    Converts a list of CMATRIX objects into a list of doubles.
    These numbers are the time-series of a specifiend component (real of imaginary)
    of a specified matrix element (i,j) of each matrix:  res[k] = H[k].get(i,j).component

    Args:
        H ( list of CMATRIX(n, m) ): time-series of the matrices
        i ( int ): row index of the matrix element of interest
        j ( int ): column index of the matrix element of interest
        component ( int ): index selecting real or imaginary component
 
            - 0: real
            - 1: imaginary

    Returns:
        list of doubles: time-series of a given matrix element's component
    
    """
    sz = len(H)
    res = []
    for k in xrange(sz):
        if component==0:
            res.append( H[k].get(i,j).real )
        elif component==1:
            res.append( H[k].get(i,j).imag )

    return res


def list2MATRIX(data):
    """Converts a list of N doubles into a MATRIX(N,1) object

    Args:
        data ( list of doubles ): data to be converted

    Returns:
        MATRIX(N,1): a matrix representation of the data

    """

    N = len(data)
    res = MATRIX(N,1)
    
    for n in xrange(N):
        res.set(n, 0, data[n])

    return res



