#*********************************************************************************
#* Copyright (C) 2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: ERGO_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for dealing with the outputs of the ErgoSCF package

.. moduleauthor:: 
       Alexey V. Akimov 
  
"""


import os
import sys
import math
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
import units



def get_mtx_matrices(filename, act_sp1=None, act_sp2=None):
    """Get the matrices printed out by the DFTB+

    Args: 
        filename ( string ): the name of the file to read. In the MatrixMarket (mtx) format
        act_sp1 ( list of ints or None): the row active space to extract from the original files
            Indices here start from 0. If set to None - the number of AOs will be
            determined automatically from the file. [default: None]
        act_sp2 ( list of ints or None): the cols active space to extract from the original files
            Indices here start from 0. If set to None - the number of AOs will be
            determined automatically from the file. [default: None]
    
    Returns: 
        list of CMATRIX(N, M): X: where N = len(act_sp1) and M = len(act_sp2) 
            These are the corresponding property matrices (converted to the complex type)
    
    """

    
    # Copy the content of the file to a temporary location
    f = open(filename, "r")
    A = f.readlines()
    f.close()    

    for a in A:
        


    norbs = int(float(A[1].split()[1]))
    nkpts = int(float(A[1].split()[2]))


    # Determine the dimensions of the output matrices
    if act_sp1==None:
        nstates1 = norbs
        act_sp1 = range(0, nstates1)
    else:
        nstates1 = len(act_sp1)

    if act_sp2==None:
        nstates2 = norbs
        act_sp2 = range(0, nstates2)
    else:
        nstates2 = len(act_sp2)

           
    # Output variables    
    X = []
    
    # For each k-point
    for ikpt in xrange(nkpts):

        # Write just the matrix data into a temporary files
        # This procedure is made fault-resistant, to avoid wrong working 
        # when some of the matrix elements printed out are nonsense.
        B = A[5+ikpt*norbs : 5+(1+ikpt)*norbs]
        
        f1 = open(filename+"_tmp", "w")  
        for i in xrange(norbs):
            
            tmp = B[i].split()
            line = ""            
            for j in xrange(norbs): 
                z = 0.0
                if tmp[j]=="NaN":
                    z = 0.0
                else:
                    try:
                        z = float(tmp[j])
                    except ValueError:
                        z = 0.0
                        pass
                    except TypeError:
                        z = 0.0
                        pass                    
                line = line + "%10.8f  " % (z)
            line = line + "\n"

            f1.write(line)
        f1.close()
                 
        # Read in the temporary file - get the entire matrix 
        x = MATRIX(norbs, norbs);  
        x.Load_Matrix_From_File(filename+"_tmp")        
        
        # Extract the sub-matrix of interest
        x_sub = MATRIX(nstates1, nstates2)
        pop_submatrix(x, x_sub, act_sp1, act_sp2)

        # Add the resulting matrices to the returned result
        X.append( CMATRIX(x_sub) )    
        
    return X




def xyz_traj2gen_sp(infile, md_iter):
    """
 
    This file converts the xyz trajectory file 
    to a string that contains the geometry of the ```md_iter``` step
    that can be used to setup the input for the ErgoSCF calculations

    Args:
        infile ( string ): the name of the input xyz trajectory file
        md_iter ( int ): index of the timeframe to extract
    
    Returns:
        string: atoms and coordinates section [ units: same as in xyz file ]
        
    """

    f = open(infile, "r")
    A = f.readlines()
    f.close()
    sz = len(A)

    # Determine the key info
    nat = int( float(A[0].split()[0] ) )

    # Make up the output file
    line = ""
    for i in xrange(nat):
        ln_indx = (nat+2)*md_iter + (i + 2)
        tmp = A[ln_indx].split()

        at = tmp[0]
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %s %15.8f  %15.8f  %15.8f\n" % (at, x, y, z)

    return line 

    
         
def xyz_traj2gen_ovlp(infile, md_iter1, md_iter2):
    """

    This file converts the xyz trajectory file 
    to a string that contains the superimposed geometries of the ```md_iter1```
    and ```md_step2``` steps  that can be used to setup the input for 
    the ErgoSCF calculations
 
    Args:
        infile ( string ): the name of the input xyz trajectory file
        md_iter1 ( int ): index of the first timeframe to extract
        md_iter2 ( int ): index of the second timeframe to extract
    
    Returns:
        string: atoms and coordinates section [ units: same as in xyz file ]
        
    """

    f = open(infile, "r")
    A = f.readlines()
    f.close()
    sz = len(A)

    # Determine the key info
    nat = int( float(A[0].split()[0] ) )

    # Make up the output file
    line = ""
    for i in xrange(nat):
        ln_indx = (nat+2)*md_iter1 + (i + 2)
        tmp = A[ln_indx].split()

        at = tmp[0]
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %s %15.8f  %15.8f  %15.8f\n" % (at, x, y, z)

    for i in xrange(nat):
        ln_indx = (nat+2)*md_iter2 + (i + 2)
        tmp = A[ln_indx].split()

        at = tmp[0]
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %s %15.8f  %15.8f  %15.8f\n" % (at, x, y, z)

    return line 

