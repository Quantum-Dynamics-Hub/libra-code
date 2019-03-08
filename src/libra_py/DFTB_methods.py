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
.. module:: DFTB_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for dealing with the outputs from DFTB+ package

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

import units
import regexlib as rgl


def get_energy_forces(filename, nat):
    """Get forces from the input file 

    Args:
        filename ( string ): the name of the file to read, usually 
            this is a the file "detailed.out" produced with the input
            containing option "CalculateForces = Yes "

        nat ( int ): the number of atoms in the system

    Returns:
        tuple: ( E_ex, F, Flst ), where:

            * E_ex ( double ): excitation energy for the given state [ units: a.u. ]
            * F ( MATRIX(ndof, 1) ): the forces acting on all atoms [ units: a.u. of force ]
            * Flst ( list of [x, y, z] ): the forces acting on all atoms [ units: a.u. of force ]

    Warning: 
        it is likely not gonna work for files other than "detailed.out"
    
    """

    # Read the file
    f = open(filename, "r")
    A = f.readlines()
    f.close()
    

    # Returned properties initialized
    E_ex = 0.0
    F = MATRIX(3 * nat, 1)
    Flst = []
    
    # Lets look for what we need
    sz = len(A)
    for i in xrange(sz):    
        tmp = A[i].split()
        
        #=========== Look for forces =============
        if len(tmp)==2:
            if tmp[0]=="Total" and tmp[1]=="Forces":
                for j in xrange(i+1, i+1+nat):
                    ind = j - i - 1 
                    tmp1 = A[j].split()
                    x = float(tmp1[0])
                    y = float(tmp1[1])
                    z = float(tmp1[2])
                    F.set(3*ind+0, 0, x); F.set(3*ind+1, 0, y); F.set(3*ind+2, 0, z)
                    Flst.append([x, y, z])
                    
        #=========== Excitation energy =============
        if len(tmp)>3:
            if tmp[0]=="Excitation" and tmp[1]=="Energy:":
                E_ex = float(tmp[2])
        
    return E_ex, F, Flst



def get_dftb_matrices(filename, act_sp1, act_sp2):
    """Get the matrices printed out by the DFTB+

    Args: 
        filename ( string ): the name of the file to read, usually 
            these are any of the files: "hamsqr1.dat", "hamsqr2.dat", etc. or 
            "oversqrt.dat". Produced with the input containing option "WriteHS = Yes"
        act_sp1 ( list of ints ): the row active space to extract from the original files
            Indices here start from 0
        act_sp2 ( list of ints ): the cols active space to extract from the original files
            Indices here start from 0
    
    Returns: 
        list of CMATRIX(N, M): X: where N = len(act_sp1) and M = len(act_sp2) 
            These are the corresponding property matrices (converted to the complex type)
            for each k-point, such that X[ikpt] is a CMATRIX(N, M) containing the 
            overlaps/Hamiltonians for the k-point ```ikpt```.

    Warning:
        So far tested only for a single k-point!
    
    """

    # Determine the dimensions of the output matrices
    nstates1 = len(act_sp1)
    nstates2 = len(act_sp2)
    
    # Copy the content of the file to a temporary location
    f = open(filename, "r")
    A = f.readlines()
    f.close()    
    norbs = int(float(A[1].split()[1]))
    nkpts = int(float(A[1].split()[2]))
           
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


def xyz_traj2gen_sp(infile, outfile, md_iter, sys_type):
    """
 
    This file converts the xyz trajectory file in the format produced
    by DFTB+ to a gen file containing the `md_iter`-th step geometry


    Args:
        infile ( string ): the name of the input xyz trajectory file
        outfile ( string ): the name of a output gen single point file
        md_iter ( int ): index of the timeframe to extract
        sys_type ( string ): "C" or "P" - the system type - cluster or periodic
    
    Returns:
        none: but creates a file
        
    """

    f = open(infile, "r")
    A = f.readlines()
    f.close()
    sz = len(A)

    # Determine the key info
    nat = int( float(A[0].split()[0] ) )
    at_types = []

    for i in xrange(nat):
        at_typ = A[i+2].split()[0]   # atom type
        if at_typ not in at_types:
            at_types.append(at_typ)


    # Make up the output file
    line = "%5i  %s \n" % (nat, sys_type)

    for typ in at_types:
        line = line + " %s " % (typ)
    line = line + "\n"

    for i in xrange(nat):
        ln_indx = (nat+2)*md_iter + (i + 2)
        tmp = A[ln_indx].split()

        at_indx = i + 1
        at = tmp[0]
        typ_indx = at_types.index(at) + 1
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %5i %3i  %15.8f  %15.8f  %15.8f\n" % (at_indx, typ_indx, x, y, z)


    # Write the file
    f1 = open(outfile, "w")
    f1.write(line)
    f1.close()

    
         
def xyz_traj2gen_ovlp(infile, outfile, md_iter1, md_iter2, sys_type):
    """
 
    This file converts the xyz trajectory file in the format produced
    by DFTB+ to a gen file containing the superimposed `md_iter1`-th 
    and `md_iter2`-th steps geometries

    Args:
        infile ( string ): the name of the input xyz trajectory file
        outfile ( string ): the name of a output gen single point file
        md_iter1 ( int ): index of the first timeframe to extract
        md_iter2 ( int ): index of the second timeframe to extract
        sys_type ( string ): "C" or "P" - the system type - cluster or periodic
    
    Returns:
        none: but creates a file
        
    """

    f = open(infile, "r")
    A = f.readlines()
    f.close()
    sz = len(A)

    # Determine the key info
    nat = int( float(A[0].split()[0] ) )
    at_types = []

    for i in xrange(nat):
        at_typ = A[i+2].split()[0]   # atom type
        if at_typ not in at_types:
            at_types.append(at_typ)


    # Make up the output file
    line = "%5i  %s \n" % (2*nat, sys_type)

    for typ in at_types:
        line = line + " %s " % (typ)
    line = line + "\n"


    for i in xrange(nat):
        ln_indx = (nat+2)*md_iter1 + (i + 2)
        tmp = A[ln_indx].split()

        at_indx = i + 1
        at = tmp[0]
        typ_indx = at_types.index(at) + 1
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %5i %3i  %15.8f  %15.8f  %15.8f\n" % (at_indx, typ_indx, x, y, z)

    for i in xrange(nat):
        ln_indx = (nat+2)*md_iter2 + (i + 2)
        tmp = A[ln_indx].split()

        at_indx = nat + i + 1
        at = tmp[0]
        typ_indx = at_types.index(at) + 1
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %5i %3i  %15.8f  %15.8f  %15.8f\n" % (at_indx, typ_indx, x, y, z)


    # Write the file
    f1 = open(outfile, "w")
    f1.write(line)
    f1.close()










