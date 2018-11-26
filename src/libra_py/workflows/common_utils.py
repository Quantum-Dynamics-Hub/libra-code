#*********************************************************************************
#* Copyright (C) 2017-2018 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  
#
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
#from libra_py import *

    
def get_matrix(nrows, ncols, filename_re, filename_im, act_sp):
    """
    This file reads the real and imaginary components of a matrix of given original size,  
    takes its sub-matrix (as defined by the act_sp function) and returns the resulting 
    complex matrix

    nrows (int)   - the number of rows in the original matrix (read from the files)
    ncols (int)   - the number of columns in the original matrix (read from the files)
    filename_re (string) - the name of the file containing the real part of the matrix 
    filename_im (string) - the name of the file containing the imaginary part of the matrix 
    act_sp (list of ints) - the indices of the columns and rows to be taken to construct the resulting matrices
    
    """

    X_re = MATRIX(nrows, ncols); X_re.Load_Matrix_From_File(filename_re)
    X_im = MATRIX(nrows, ncols); X_im.Load_Matrix_From_File(filename_im)

    nstates = len(act_sp)
    x_re = MATRIX(nstates, nstates);
    x_im = MATRIX(nstates, nstates);

    pop_submatrix(X_re, x_re, act_sp, act_sp)
    pop_submatrix(X_im, x_im, act_sp, act_sp)

    return CMATRIX(x_re, x_im)

