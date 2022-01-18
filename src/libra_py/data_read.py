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
import scipy.sparse as sp

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import common_utils as comn
from libra_py import data_conv
import util.libutil as comn

    
def get_matrix(nrows, ncols, filename_re, filename_im, act_sp, get_real=1, get_imag=1):
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
        get_real ( int ): whether we want to read the real component [ default: 1 - Yes ]
        get_imag ( int ): whether we want to read the imaginary component [ default: 1 - Yes ]

    Returns:
        CMATRIX(N, N): where N is the number of actively included rows/columns


    Example:

        The following snippet will create a 4 x 4 matrix from the files "Ham_0_re"
        and "Ham_0_im" from the "res" directory. Each of the files is expected to be 
        a matrix of 10 x 10 in size. The rezulting 4 x 4 matrix will contain entries on
        the intersection of columns and rows with indices 0, 1, 3, and 4.

        >>> X = get_matrix(10, 10, "res/Ham_0_re", "res/Ham_0_im", [0,1,3,4])


    """
   
    X_re = MATRIX(nrows, ncols); 
    if get_real==1:
        if os.path.exists(filename_re):
            X_re.Load_Matrix_From_File(filename_re)
        else:
            print(F"File {filename_re} does not exist. Exiting...")
            sys.exit(0)
   
    X_im = MATRIX(nrows, ncols); 
    if get_imag==1:
        if os.path.exists(filename_re):       
            X_im.Load_Matrix_From_File(filename_im)
        else:
            print(F"File {filename_im} does not exist. Exiting...")
            sys.exit(0)


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
            * **params["get_real"]** ( int ): whether we want to read the real component [ default: 1 - Yes ]
            * **params["get_imag"]** ( int ): whether we want to read the imaginary component [ default: 1 - Yes ]
            

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
    default_params = { "data_re_suffix":"_re", "data_im_suffix":"_im", "active_space":range(params["data_dim"]), "get_real":1, "get_imag":1}
    comn.check_input(params, default_params, critical_params)

    ndim = params["data_dim"]  # the number of cols/row in the input files

    data = []
    for i in range(params["isnap"],params["fsnap"]):

        filename_re = params["data_re_prefix"]+str(i)+params["data_re_suffix"]
        filename_im = params["data_im_prefix"]+str(i)+params["data_im_suffix"]
        data_i = get_matrix(ndim, ndim, filename_re, filename_im, params["active_space"], params["get_real"], params["get_imag"]) 
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
                In addition, requires parameters described in :func:`get_data`

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

   
def read_2D_grid(filename):
    """
    This function reads the 2D map pyplot formatted data from a file.
    
    Args:
        filename (sting) : name of the file to read
        
    Returns:
        double list, double list, list of lists of doubles:
            X grid, Y grid, and the Z values of the grid points
    
    """
    x,y,z = [], [], []
    
    f = open(filename, "r")
    A = f.readlines()
    f.close()
    
    #========== Determine the numbers =====
    nlines = len(A)        
    ny = 0
    line_size = len(A[ny].split())    
    while line_size > 0:        
        ny += 1
        line_size = len(A[ny].split())
                        
    nx = int((nlines + 1)/(ny + 1))
        
    #========== Get the grids =====
    for ix in range(nx):
        x.append( float(A[ix*(ny+1)].split()[0]) )
        
    for iy in range(ny):
        y.append( float(A[iy].split()[1]) )
    
    for ix in range(nx):
        z_x = []
        for iy in range(ny):
            z_xy = float(A[ix*(ny+1)+iy].split()[2])
            z_x.append(z_xy)
        z.append(z_x)
        
    return x, y, z
          
  
def get_all_data(params):
    """
    This function reads all the data needed for running the nonadiabatic dynamics.
    Args:
        params (dictionary):
            * **params['istep']*: The initial step
            * **params['nsteps']*: Number of steps to be read from the inital step
            * **params['path_to_res_files']*: The full path to where the data are stored
            * **params['Hvib_re_prefix']*: The prefix of the real part of the Hamiltonian file
            * **params['Hvib_im_prefix']*: The prefix of the imaginary part of the Hamiltonian file
            * **params['Hvib_re_suffix']*: The suffix of the real part of the Hamiltonian file
            * **params['Hvib_im_suffix']*: The suffix of the imaginary part of the Hamiltonian file
            * **params['St_re_prefix']*: The prefix of the real part of the time-overlap file
            * **params['St_re_suffix']*: The suffix of the real part of the time-overlap file
            * **params['is_sparse']**: whether the files are stored as sparse (`npz`) or readable text files
    Returns:
        Hadi (list of CMATRIX): A list of real part of the Hamiltonian
        Hvib (list of CMATRIX): A list of the full Hamiltonian
        nac (list of CMATRIX): A list of the nonadiabatic couplings
        St (list of CMATRIX): A list of the time-overlap matrices
        nstates (integer): The number of states which is equivalent to the number of rows of one of the files.
    """
    Hadi, Hvib, nac, St = [], [], [], []
    istep = params['istep']
    fstep = istep + params['nsteps']
    path_to_res_files = params['path_to_res_files']
    hvib_re_prefix = params['Hvib_re_prefix']
    hvib_im_prefix = params['Hvib_im_prefix']
    st_re_prefix   = params['St_re_prefix']
    hvib_re_suffix = params['Hvib_re_suffix']
    hvib_im_suffix = params['Hvib_im_suffix']
    st_re_suffix   = params['St_re_suffix']
    # find nstates
    if params['is_sparse']:
        tmp_mat = sp.load_npz(F'{path_to_res_files}/{hvib_re_prefix}{istep}{hvib_re_suffix}.npz')
        nstates = tmp_mat.shape[0]
    else:
        tmp_mat = np.loadtxt(F'{path_to_res_files}/{hvib_re_prefix}{istep}{hvib_re_suffix}')
        nstates = tmp_mat.shape[0]
    dummy = MATRIX(nstates, nstates)
    if params['is_sparse']:
        for step in range(istep, fstep):
            print('Reading the data of step', step)
            filename = F'{path_to_res_files}/{hvib_re_prefix}{step}{hvib_re_suffix}.npz'
            hvib_re = np.array(sp.load_npz(filename).todense().real)
            hvib_re = data_conv.nparray2MATRIX(hvib_re)
            filename = F'{path_to_res_files}/{hvib_im_prefix}{step}{hvib_im_suffix}.npz'
            hvib_im = np.array(sp.load_npz(filename).todense().real)
            hvib_im = data_conv.nparray2MATRIX(hvib_im)
            filename = F'{path_to_res_files}/{st_re_prefix}{step}{st_re_suffix}.npz'
            st_re = np.array(sp.load_npz(filename).todense().real)
            st_re = data_conv.nparray2MATRIX(st_re)
            Hvib.append( CMATRIX(hvib_re, hvib_im) )
            Hadi.append( CMATRIX(hvib_re, dummy) )
            nac.append( CMATRIX(-1.0*hvib_im, dummy) )
            St.append( CMATRIX(st_re, dummy) )
    else:
        for step in range(istep, fstep):
            # Hvib
            print('Reading the data of step', step)
            filename_re = F'{path_to_res_files}/{hvib_re_prefix}{step}{hvib_re_suffix}'
            filename_im = F'{path_to_res_files}/{hvib_im_prefix}{step}{hvib_im_suffix}'
            hvib = data_read.get_matrix(nstates, nstates, filename_re, filename_im,  list(range(nstates)), 1, 1)
            Hvib.append(hvib)

            # Hadi
            Hadi.append( CMATRIX(hvib.real(), dummy) )

            # NAC
            nac.append( CMATRIX(-1.0*hvib.imag(), dummy) )

            # St
            filename_re = F'{path_to_res_files}/{st_re_prefix}{step}{st_re_suffix}'
            filename_im = F'{path_to_res_files}/{st_im_prefix}{step}{st_im_suffix}'
            st = data_read.get_matrix(nstates, nstates, filename_re, filename_im, list(range(nstates)), 1, 1)
            St.append(st)

    return Hadi, Hvib, nac, St, nstates

 
