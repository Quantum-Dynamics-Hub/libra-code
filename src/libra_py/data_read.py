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
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import common_utils as comn
import util.libutil as comn

    
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

    pop_submatrix(X_re, x_re, list(act_sp), list(act_sp))
    pop_submatrix(X_im, x_im, list(act_sp), list(act_sp))

    return CMATRIX(x_re, x_im)




def get_data(params):
    """Read a single set of data files 

    Args:
        params ( dictionary ): parameters controlling the function execution

            Required parameter keys:

            * **params["data_dim"]** ( int ): matrix dimension how many lines/columns in the file [Required!]
            * **params["active_space"]** ( list of ints ): the indices of the states we care 
                about. These indices will be used to determine the size of the created CMATRIX objects
                and only these states will be extracted from the original files [ default: range(data_dim) ]
            * **params["isnap"]** ( int ): index of the first file to read [Required!]
            * **params["fsnap"]** ( int ): index of the final file to read [Required!]
            * **params["data_re_prefix"]** ( string ): prefixes of the files with real part of the data [Required!]
            * **params["data_im_prefix"]** ( string ): prefixes of the files with imaginary part of the data [Required!]
            * **params["data_re_suffix"]** ( string ): suffixes of the files with real part of the Hvib(t) [default: "_re"]
            * **params["data_im_suffix"]** ( string ): suffixes of the files with imaginary part of the Hvib(t) [default: "_im"]

    Returns:
        list of CMATRIX objects: data: 
            a time series of data matrices, such that data[time] is a data at time step `time`

    Example:
        This example will read 10 pairs of files: "Hvib_0_re", "Hvib_0_im", "Hvib_1_re", "Hvib_1_im", ...
        "Hvib_9_re", "Hvib_9_im". Each file should contain a 4 x 4 matrix of numbers. It will generate a 
        list of 4 x 4 complex-valued matrices.

        >>> hvib = get_data({"data_dim":4, "isnap":0, "fsnap":10, "data_re_prefix":"Hvib", "data_im_prefix":"Hvib"})

        The following example will do the same as the example above, however the intially-read 4 x 4 matrices will
        be partially discarded. Out of 16 values only 4 (the upper left block of 4 numbers)  will be stored in 
        the resulting list of 2 x 2 complex-valued matrices. 

        >>> hvib = get_data({"data_dim":4, "isnap":0, "fsnap":10, "data_re_prefix":"Hvib", "data_im_prefix":"Hvib", "active_space":[0,1]})


    """

    critical_params = ["data_dim", "isnap", "fsnap", "data_re_prefix", "data_im_prefix"]
    default_params = { "data_re_suffix":"_re", "data_im_suffix":"_im", "active_space":range(params["data_dim"])}
    comn.check_input(params, default_params, critical_params)

    ndim = params["data_dim"]  # the number of cols/row in the input files

    data = []
    for i in range(params["isnap"],params["fsnap"]):

        filename_re = params["data_re_prefix"]+str(i)+params["data_re_suffix"]
        filename_im = params["data_im_prefix"]+str(i)+params["data_im_suffix"]
        data_i = get_matrix(ndim, ndim, filename_re, filename_im, params["active_space"] ) 
        data.append(data_i)

    return data


def get_data_sets(params):
    """Reads several sets of data files 

    Args:
        params ( dictionary ): parameters controlling the function execution [Required!]

            Required parameter keys:

            * **params["data_set_paths"]** ( list of strings ):
                define the paths of the directories where the data files for
                different data sets (e.g. independent MD trajectories) are located. 
            .. note::
                In addition, requires parameters described in
                :func:`get_data`

    Returns:
        list of lists of CMATRIX: data: 
            the time series of Hvib matrices for several data sets, such that
            data[idata][time] is a CMATRIX for the data set indexed by `idata`
            at time `time`


    Example:
        The full name of the vibronic Hamiltonian files read by this module should be:
    
        params["data_set_paths"][idata]+params["data_re_prefix"]+integer(time step)+params["data_re_suffix"] - for real part

        params["data_set_paths"][idata]+params["data_im_prefix"]+integer(time step)+params["data_im_suffix"] - for imaginary part

        Say, the directory "/home/alexeyak/test/step3/res0" contains files:
        Hvib_0_re, Hvib_1_re, .... ,    Hvib_999_re
        Hvib_0_im, Hvib_1_im, .... ,    Hvib_999_im

        Then set:

        >>> params["data_set_paths"] = ["/home/alexeyak/test/step3/res0/"]
        >>> params["data_re_prefix"] = "Hvib_"
        >>> params["data_re_suffix"] = "_re"
        >>> params["data_im_prefix"] = "Hvib_"
        >>> params["data_im_suffix"] = "_im"

    """

    critical_params = [ "data_set_paths" ] 
    default_params = { }
    comn.check_input(params, default_params, critical_params)

    data = []

    for idata in params["data_set_paths"]:   # over all MD trajectories (data sets)
        prms = dict(params)    
        prms.update({"data_re_prefix": idata+params["data_re_prefix"] })
        prms.update({"data_im_prefix": idata+params["data_im_prefix"] })                

        data_i = get_data(prms)  
        data.append(data_i)

    return data






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




def get_data_from_file2(filename, cols):
    """Read in the numeric data stored in a file as columns into Python lists

    Args:
        filename ( string ): The name of the data file
        cols ( list of ints ): the indices of the columns to read

    Returns:
        (list of lists): data

    """

    f = open(filename,"r")
    A = f.readlines()
    f.close()

    sz = len(cols)
    res = []
    for i in range(0,sz):
        res.append([])

    for a in A:
        tmp = a.split()

        for i in range(0,sz):

            x = float(tmp[ cols[i] ])
            res[i].append(x)

    return res

   
  
   
