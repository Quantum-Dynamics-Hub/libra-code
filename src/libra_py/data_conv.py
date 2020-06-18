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
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import common_utils as comn
import util.libutil as comn


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
        None: but transforms X input directly, so changes the original input


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

    for i in range(0,nstates):
        for j in range(0,nstates):
            scl.set(i,j, 1.0+0.0j)

    critical_params = [  ] 
    default_params = { "shift1":sh1, "shift2":sh2, "scale":scl  }
    comn.check_input(params, default_params, critical_params)

    for idata in range(0,ndata):
        for istep in range(0,nsteps):

            tmp = CMATRIX(X[idata][istep])
            tmp = tmp + params["shift1"]
            tmp.dot_product( tmp, params["scale"] )
            tmp = tmp + params["shift2"]
            X[idata][istep] = CMATRIX(tmp)



def unit_conversion(X, scaling_factor):
    """

    Rescales the data X uniformly:     
    X(original) ->    ( X ) (x) (scaling_factor),  

    Here, (x) indicates the element-wise multiplicaiton, and shift1, shift2, and scale are matrices

    Args: 
        X ( list of lists of CMATRIX(nstates, nstates) ): the original data stored as
            X[idata][step] - a CMATRIX(nstates, nstates) for the dataset `idata` and time 
            step `step`
        scaling_factor ( double or complex ): the rescaling factor

    Returns:    
        None: but transforms X input directly, so changes the original input

    """

    nst = X[0][0].num_of_cols

    scl = CMATRIX(nst, nst);                             
    shi = CMATRIX(nst, nst);                             

    # Default ones
    for i in range(0,nst):
        for j in range(0,nst):
            scl.set(i,j, scaling_factor*(1.0+0.0j))
            shi.set(i,j, (0.0+0.0j))        
    transform_data(X, {"shift2":shi, "scale":scl }) 




def scale_NAC(X, a, b, scaling_factor):
    """

    Rescales only the matrix elements X_ab by a uniform scaling factor:     
    X_ab(original) ->    ( X_ab ) x (scaling_factor),  

    Args: 
        X ( list of lists of CMATRIX(nstates, nstates) ): the original data stored as
            X[idata][step] - a CMATRIX(nstates, nstates) for the dataset `idata` and time 
            step `step`
        scaling_factor ( double or complex ): the rescaling factor

    Returns:    
        None: but transforms X input directly, so changes the original input

    """


    nst = X[0][0].num_of_cols

    scl = CMATRIX(nst, nst);                             
    shi = CMATRIX(nst, nst);                             

    # Default ones
    for i in range(0,nst):
        for j in range(0,nst):
            scl.set(i,j, 1.0+0.0j)
            shi.set(i,j, 0.0+0.0j)        
    scl.set(a,b, scaling_factor)

    transform_data(X, {"shift2":shi, "scale":scl }) 



def scale_NACs(X, scaling_factor):
    """

    Rescales all the off-diagonal matrix elements X_ab by a uniform scaling factor:     
    X_ab(original) ->    ( X_ab ) x (scaling_factor),  for all a!=b

    Args: 
        X ( list of lists of CMATRIX(nstates, nstates) ): the original data stored as
            X[idata][step] - a CMATRIX(nstates, nstates) for the dataset `idata` and time 
            step `step`
        scaling_factor ( double or complex ): the rescaling factor

    Returns:    
        None: but transforms X input directly, so changes the original input

    """


    nst = X[0][0].num_of_cols

    scl = CMATRIX(nst, nst);                             
    shi = CMATRIX(nst, nst);                             

    # Default ones
    for i in range(0,nst):
        for j in range(0,nst):
            if i!=j:
                scl.set(i,j, scaling_factor*(1.0+0.0j))
            else:
                scl.set(i,j, (1.0+0.0j))
            shi.set(i,j, (0.0+0.0j))        

    transform_data(X, {"shift2":shi, "scale":scl }) 




def scissor(X, a, dE):
    """

    Shift the diagonal elements of X : X_ii for all i = a, a+1, ... by a constant value dE
    X_ii(original) ->    ( X_ii ) + dE, for all i >= a

    Args: 
        X ( list of lists of CMATRIX(nstates, nstates) ): the original data stored as
            X[idata][step] - a CMATRIX(nstates, nstates) for the dataset `idata` and time 
            step `step`
        dE ( double or complex ): the shift magnitude

    Returns:    
        None: but transforms X input directly, so changes the original input

    """


    nst = X[0][0].num_of_cols

    scl = CMATRIX(nst, nst);                             
    shi = CMATRIX(nst, nst);                             

    # Default ones
    for i in range(0,nst):
        for j in range(0,nst):
            scl.set(i,j, 1.0+0.0j)
            shi.set(i,j, 0.0+0.0j)        

    for i in range(a, nst):
        shi.add(i,i, dE)

    transform_data(X, {"shift2":shi, "scale":scl }) 




def unpack1(H, i, j, component=2):
    """

    Converts a list of CMATRIX or MATRIX objects into a list of doubles.
    These numbers are the time-series of a specifiend component (real of imaginary)
    of a specified matrix element (i,j) of each matrix:  res[k] = H[k].get(i,j).component

    Args:
        H ( list of CMATRIX(n, m) or MATRIX(n, m)): time-series of the matrices
        i ( int ): row index of the matrix element of interest
        j ( int ): column index of the matrix element of interest
        component ( int ): index selecting real or imaginary component
 
            - 0: real
            - 1: imaginary
            - 2: the whole thing - use for real matrices [ default ]

    Returns:
        list of doubles: time-series of a given matrix element's component
    
    """
    sz = len(H)
    res = []
    for k in range(0,sz):
        if component==0:
            res.append( H[k].get(i,j).real )
        elif component==1:
            res.append( H[k].get(i,j).imag )
        elif component==2:
            res.append( H[k].get(i,j) )

    return res



def unpack2(H, i, component=2):
    """

    Converts a CMATRIX or MATRIX object into a list of doubles.
    These numbers are the time-series of a specifiend component (real of imaginary)
    of a specified matrix column ```i``` res[k] = H.get(k,i).component

    Args:
        H ( CMATRIX(n, m) or MATRIX(n, m) ): time-series in a form of a matrix
        i ( int ): column index of interest
        component ( int ): index selecting real or imaginary component
 
            - 0: real
            - 1: imaginary
            - 2: the whole thing - use for real matrices

    Returns:
        list of doubles: time-series of a given matrix column's components
    
    """

    sz = H.num_of_rows

    res = []
    for k in range(0,sz):
        if component==0:
            res.append( H.get(k,i).real )
        elif component==1:
            res.append( H.get(k,i).imag )
        elif component==2:
            res.append( H.get(k,i) )

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
    
    for n in range(0,N):
        res.set(n, 0, data[n])

    return res


def nparray2MATRIX(data):
    """
    Converts 2D np.array of shape( N, M ) doubles into a MATRIX( N, M ) object
    The numpy array can contain either complex or real values    

    Args:
        data ( 2D np.array of dimension N x M ): data to be converted
    Returns:
        MATRIX( N, M ): a matrix representation of the data
    
    """

    N = data.shape[0]
    M = data.shape[1]
    
    res = MATRIX(N,M)
    
    for n in range(0,N):
        for m in range(0,M):
            res.set(n, m, data[n][m].real)

    return res


def nparray2CMATRIX(data):
    """
    Converts 2D np.array of shape( N, M ) doubles into a CMATRIX( N, M ) object
    The numpy array should be complex   
 
    Args:
        data ( 2D np.array of dimension N x M ): data to be converted
    Returns:
        MATRIX( N, M ): a matrix representation of the data
    
    """

    N = data.shape[0]
    M = data.shape[1]

    res = CMATRIX(N,M)

    for n in range(0,N):
        for m in range(0,M):
            res.set(n, m, data[n][m])

    return res



def matrix2list(q):
    """
    Converts the MATRIX(ndof, 1) or CMATRIX(ndof, 1) to a list of `ndof` float/complex numbers

    Args:

=======
        q ( MATRIX(ndof, 1) or CMATRIX(ndof, 1) ): input matrix


    Returns:
   
        list : list representation of the matrix

    """
    
    list_q = []
    
    ndof = q.num_of_rows
    
    for idof in range(ndof):
        list_q.append( q.get(idof, 0) )
        
    return list_q



def make_list(nitems, value):
    """
    Creates a list of `nitems` items, each of which is `value`

    Args:
        nitems ( int ): the size of the resulting list to create
        value ( any type ): the value of each initialized element of the list
    
    Returns:
        list : the list of the added values

    """
    
    res_list = []
    
    for iitem in range(nitems):
        res_list.append( value )
        
    return res_list

  
  
def form_block_matrix( mat_a, mat_b, mat_c, mat_d ):
    """
    This function gets four numpy arrays and concatenate them into a 
    new matrix in a block format. These matrices should have the same 
    shape on each side which they get concatenated.

         |mat_a   mat_b|
    S =  |             |
         |mat_c   mat_d|

    Args:

        mat_a, mat_b, mat_c, mat_d ( 2D numpy arrays ): The matrices which will form the block matrix

    Returns:

        block_matrix ( numpy 2D array ): The block matrix of the four matrices above.

    """

    # Concatenate the two marix on their row axis
    block_1 = np.concatenate( ( mat_a, mat_b ) )
    block_2 = np.concatenate( (mat_c, mat_d ) )

    # Now concatenate the above concatenated matrices on their 
    # column axis and form the final matrix
    block_matrix = np.concatenate( ( block_1, block_2 ), axis=1 )

    return block_matrix

