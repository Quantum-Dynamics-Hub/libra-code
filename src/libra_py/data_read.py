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
.. module:: data_read
   :platform: Unix, Windows
   :synopsis: 
       This module implements various functions for getting data from files

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

    
def get_matrix(nrows, ncols, filename_re, filename_im, act_sp):
    """

    This file reads the real and imaginary components of a matrix of given original size,  
    takes its sub-matrix (as defined by the act_sp function) and returns the resulting 
    complex matrix

    Args:
        nrows ( int ): the number of rows in the original matrix (read from the files)
        ncols ( int ): the number of columns in the original matrix (read from the files)
        filename_re ( string ): the name of the file containing the real part of the matrix 
        filename_im ( string ): the name of the file containing the imaginary part of the matrix 
        act_sp ( list of N ints ): the indices of the columns and rows to be taken to construct the 
            resulting matrices. The indexing starts from 0. These numbers shold not be larger than 
            `nrows` or `ncols`

    Returns:
        CMATRIX(N, N): where N is the number of actively included rows/columns


    Example:

        The following snippet will create a 4 x 4 matrix from the files "Ham_0_re"
        and "Ham_0_im" from the "res" directory. Each of the files is expected to be 
        a matrix of 10 x 10 in size. The rezulting 4 x 4 matrix will contain entries on
        the intersection of columns and rows with indices 0, 1, 3, and 4.

        >>> X = get_matrix(10, 10, "res/Ham_0_re", "res/Ham_0_im", [0,1,3,4])


    """

    X_re = MATRIX(nrows, ncols); X_re.Load_Matrix_From_File(filename_re)
    X_im = MATRIX(nrows, ncols); X_im.Load_Matrix_From_File(filename_im)

    nstates = len(act_sp)
    x_re = MATRIX(nstates, nstates);
    x_im = MATRIX(nstates, nstates);

    pop_submatrix(X_re, x_re, act_sp, act_sp)
    pop_submatrix(X_im, x_im, act_sp, act_sp)

    return CMATRIX(x_re, x_im)




def get_data_from_file(filename, xindx, yindx, xminval=None, xmaxval=None, yminval=None, ymaxval=None):
    """Read in the numeric data stored in a file as columns into Python lists

    Args:
        filename ( string ): The name of the data file
        xindx ( int ): the index of the column read as X
        yindx ( int ): the index of the column read as Y
        xminval ( double ): the minimal X value allowed in the read data set, 
            the points with X values below it will not be included [ default: None ]
        xmaxval ( double ): the maximal X value allowed in the read data set, 
            the points with X values above it will not be included [ default: None ]
        yminval ( double ): the minimal Y value allowed in the read data set, 
            the points with Y values below it will not be included [ default: None ]
        ymaxval ( double ): the maximal Y value allowed in the read data set, 
            the points with Y values above it will not be included [ default: None ]

    Returns:
        (list, list): (X, Y), where: 

            * X ( list of doubles ): x values read from the file, cropped according the conditions
            * Y ( list of doubles ): y values read from the file, cropped according the conditions

    """

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


