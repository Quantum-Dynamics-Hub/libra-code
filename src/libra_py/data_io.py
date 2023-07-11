#*********************************************************************************                     
#* Copyright (C) 2019-2020 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: data_io
   :platform: Unix, Windows
   :synopsis: This module implements the read/write functions specifically designed to work with some data types
       List of functions:
           * add_intlist2file(filename, t, X)
           * add_doublelist2file(filename, t, X)
           * add_matrix2file(filename, t, X)
           * add_cmatrix2file(filename, t, X)
           * file2intlist(filename)
           * file2doublelist(filename)
           * file2matrix(filename, nrows, ncols)
           * file2matrix(filename, nrows, ncols)

.. moduleauthor:: Alexey V. Akimov
  
"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"


import os
import sys
import math
import copy
import re
import numpy as np
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
from . import units
from . import data_read



def add_intlist2file(filename, t, X):
    """
    This function appends a new line of type: [t, X[0], X[1], ... X[sz-1] ] to a file. 
    Where sz = len(X)

    Args:
        filename ( string ): the name of the file to where the data will be printed out
        t ( float ): usually the time marker for the output
        X ( list of ints ): object containing the data to be printed out

    Returns:
        None: but updates the existing file

    """
            
    f = open(filename, "a")
    
    line = "%5.3f" % (t)
    for x in X:        
        line = line + " %5i" % (x)
    line = line + "\n"
    
    f.write(line)
    f.close()

    
def add_doublelist2file(filename, t, X):
    """
    This function appends a new line of type: [t, X[0], X[1], ... X[sz-1] ] to a file. 
    Where sz = len(X)

    Args:
        filename ( string ): the name of the file to where the data will be printed out
        t ( float ): usually the time marker for the output
        X ( list of floats ): object containing the data to be printed out

    Returns:
        None: but updates the existing file

    """
            
    f = open(filename, "a")
    
    line = "%5.3f" % (t)
    for x in X:        
        line = line + " %8.5f" % (x)
    line = line + "\n"
    
    f.write(line)
    f.close()    


    
def add_matrix2file(filename, t, X):
    """
    This function appends a new line of type: 
    [t, X(0,0), X(0, 1), ... , X(0, ncols-1), X(1, 0), X(1, 1), ..., X(1, ncols-1), ... ] to a file. 
    Where ncols - the number of columns of the matrix X

    Args:
        filename ( string ): the name of the file to where the data will be printed out
        t ( float ): usually the time marker for the output
        X ( MATRIX(N, M) ): object containing the data to be printed out

    Returns:
        None: but updates the existing file

    """

    
    nrows = X.num_of_rows
    ncols = X.num_of_cols
    
    f = open(filename, "a")
    
    line = "%5.3f" % (t)
    for a in range(nrows):
        for b in range(ncols):
            line = line + " %8.5f" % (X.get(a,b))
    line = line + "\n"
    
    f.write(line)
    f.close()
    


def add_cmatrix2file(filename, t, X):
    """
    This function appends a new line of type: 
    [t, X(0,0).real, X(0,0).imag, X(0, 1).real, X(0, 1).imag, ... , X(0, ncols-1).real, X(0, ncols-1).imag,
        X(1,0).real, X(1,0).imag, X(1, 1).real, X(1, 1).imag, ... , X(1, ncols-1).real, X(1, ncols-1).imag,
      ... ] to a file. 
    Where ncols - the number of columns of the matrix X

    Args:
        filename ( string ): the name of the file to where the data will be printed out
        t ( float ): usually the time marker for the output
        X ( CMATRIX(N, M) ): object containing the data to be printed out

    Returns:
        None: but updates the existing file

    """

    
    nrows = X.num_of_rows
    ncols = X.num_of_cols
    
    f = open(filename, "a")
    
    line = "%5.3f" % (t)
    for a in range(nrows):
        for b in range(ncols):
            line = line + " %8.5f %8.5f" % (X.get(a,b).real, X.get(a,b).imag)
    line = line + "\n"
    
    f.write(line)
    f.close()




def file2intlist(filename):
    """
    This function reads a number of lines of type: [t, X[0], X[1], ... X[sz-1] ] from a file. 
    Where sz = len(X)
    and stores it as a list of lists of ints

    Reading analog of the :func:`libra_py.dynamics_io.add_intlist2file`

    Args:
        filename ( string ): the name of the file to where the data will be printed out

    Returns:
        list of lists of ints : object containing the data to be printed out

    """
            
    f = open(filename, "r")
    A = f.readlines()
    f.close()
    
    res = []
    for a in A:
        tmp = a.split()
        sz = len(tmp)

        res_i = []
        for i in range(1,sz):
            res_i.append( int(float(tmp[i])) )
        
        res.append(res_i)

    return res




def file2doublelist(filename):
    """
    This function reads a number of lines of type: [t, X[0], X[1], ... X[sz-1] ] from a file. 
    Where sz = len(X)
    and stores it as a list of lists of ints

    Reading analog of the :func:`libra_py.dynamics_io.add_doublelist2file`

    Args:
        filename ( string ): the name of the file to where the data will be printed out

    Returns:
        list of lists of doubles : object containing the data to be printed out

    """
            
    f = open(filename, "r")
    A = f.readlines()
    f.close()
    
    res = []
    for a in A:
        tmp = a.split()
        sz = len(tmp)

        res_i = []
        for i in range(1,sz):
            res_i.append( float(tmp[i]) )
        
        res.append(res_i)

    return res



def file2matrix(filename, nrows, ncols):
    """
    This function reads a number of lines of type: 
    [t, X(0,0), X(0, 1), ... , X(0, ncols-1), X(1, 0), X(1, 1), ..., X(1, ncols-1), ... ] from a file. 
    Where ncols - the number of columns of the matrix X

    Reading analog of the :func:`libra_py.dynamics_io.add_matrix2file`

    Args:
        filename ( string ): the name of the file to where the data will be printed out
        nrows ( int ) : the number of rows in each matrix that is stored in this file (at every time step!)
        ncols ( int ) : the number of columns in each matrix that is stored in this file (at every time step!)

    Returns:
        list of MATRIX objects : object containing the data to be printed out

    """
            
    f = open(filename, "r")
    A = f.readlines()
    f.close()    
    
    res = []
    for item in A:
        tmp = item.split()
        sz = len(tmp)

        if ncols * nrows != sz-1:
            print("ERROR in file2matrix: the dimentions of the matrix to be read do not agree with your expectations\n")
            print("The number of columns per line should be %5i, but we've got %5i" % (ncols * nrows + 1, sz) )
            print("ncols = ", ncols)
            print("nrows = ", nrows)
            print(item)
            print(tmp)
            sys.exit(0)

        res_i = MATRIX(nrows, ncols)
        cnt = 1
        for a in range(nrows):
            for b in range(ncols):
                res_i.set(a, b, float(tmp[cnt]))
                cnt += 1        
        res.append(res_i)

    return res



def file2cmatrix(filename, nrows, ncols):
    """
    This function reads a number of lines of type: 
    [t, X(0,0).real, X(0,0).imag, X(0, 1).real, X(0, 1).imag, ... , X(0, ncols-1).real, X(0, ncols-1).imag,
        X(1,0).real, X(1,0).imag, X(1, 1).real, X(1, 1).imag, ... , X(1, ncols-1).real, X(1, ncols-1).imag,
      ... ] from a file. 

    Where ncols - the number of columns of the matrix X

    Reading analog of the :func:`libra_py.dynamics_io.add_cmatrix2file`

    Args:
        filename ( string ): the name of the file to where the data will be printed out
        nrows ( int ) : the number of rows in each matrix that is stored in this file (at every time step!)
        ncols ( int ) : the number of columns in each matrix that is stored in this file (at every time step!)

    Returns:
        list of CMATRIX objects : object containing the data to be printed out

    """
            
    f = open(filename, "r")
    A = f.readlines()
    f.close()    
    
    res = []
    for item in A:
        tmp = item.split()
        sz = len(tmp)

        if 2 * ncols * nrows != sz-1:
            print("ERROR in file2cmatrix: the dimentions of the matrix to be read do not agree with your expectations\n")
            print("The number of columns per line should be %5i, but we've got %5i" % (2 * ncols * nrows + 1, sz) )
            print("ncols = ", ncols)
            print("nrows = ", nrows)
            print(item)
            print(tmp)

            sys.exit(0)

        res_i = CMATRIX(nrows, ncols)
        cnt = 1
        for a in range(nrows):
            for b in range(ncols):
                re = float(tmp[cnt])
                im = float(tmp[cnt+1])
                res_i.set(a, b,  re + 1j*im)
                cnt += 2

        res.append(res_i)

    return res


def replace_pattern_in_file(filename, old_string, new_string, method=0):

    """
    This function replaces a string with a new strin. It is mostly used to replace the 
    digits of zeros in a file to reduce the file size and its efficient storage. It is 
    the equivalent of 'sed' command in Linux.

    Args:

        filename (string): The file name.

        old_string (string): Old pattern to be replaced.

        new_string (bool): New pattern to be replaced. 

        method (int): the selector of the algorithm to do the replacement [ default: 0]
            0 - read all the file into intermediate lines
            1 - read the file without making intermediate lists, must be more memory-efficient
 

    Returns:

        None
    """

    if method==0:

        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        f = open(filename,'w')
        for i in range(len(lines)):
            f.write(re.sub(old_string, new_string, lines[i]))
        f.close()


    elif method==1:

        pattern = re.compile(re.escape(old_string), 0)  

        f = open(filename, 'r+')
        content = f.read()
        content = pattern.sub(new_string, content)
        f.seek(0)
        f.truncate()
        f.write(content)
        f.close()


def sort_hvib_file_names(files):
    """
    This function sorts the file names obtained from `glob` for the files
    Hvib_{basis}_{step}_im.npz, Hvib_{basis}_{step}_re.npz, or St_{basis}_{step}_re.npz. 
    Args:
        files (list): A list of files obtained form `glob.glob` function.
    Returns:
        files_sorted (list): The list of files which are sorted based on the step.
    """
    steps = []
    for file in files:
        digits = []
        for i in range(len(file)):
            try:
                int(file[i])
                digits.append(file[i])
            except:
                pass
        number = ''
        for i in range(len(digits)):
            number += digits[i]

        number = int(number)
        steps.append(number)
    steps = np.argsort(np.array(steps))
    files_sorted = []
    for i in range(len(steps)):
        files_sorted.append(files[steps[i]])
    return files_sorted

