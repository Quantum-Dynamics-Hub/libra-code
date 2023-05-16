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
from libra_py import units
from libra_py import data_read


def find_last_file(prefix, suffix):
    """

    This functions searches for the last file in a series of files named
    `prefix`%i`suffix`, for instance "A1.txt", "A2.txt", "A3.txt", etc.

    Args:
        prefix ( string ): the prefix of the filenames 
        suffix ( string ): the suffix of the filenames 

    Return:
        (int, string): (i, filename), where:

            * i : the largest index of the existing  file
            * filename : the name of the last existing file

    """
    stop = False
    i = 1
    while stop==False:
        exists = os.path.isfile('%s%i%s' % (prefix, i, suffix))
        if exists:
            i = i + 1
        else:
            stop = True
            
    return i-1, '%s%i%s' % (prefix, i-1, suffix)


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

    sz = len(A)       
    start = 0
    for i in range(0,sz):
        if A[i][0]!="%":
            start = i        
            break
        
    # Initialize the matrix dimentions
    tmp = A[start].split()    
    N = int(float(tmp[0]))
    M = int(float(tmp[1]))
        
    
    # Determine the dimensions of the output matrices
    nstates1, nstates2 = -1, -1
    
    if act_sp1==None:
        nstates1 = N
        act_sp1 = list(range(0, nstates1))
    else:
        nstates1 = len(act_sp1)

    if act_sp2==None:
        nstates2 = N
        act_sp2 = list(range(0, nstates2))
    else:
        nstates2 = len(act_sp2)
        
    X = CMATRIX(N,M)    
        
    for n in range(start+1, sz):
        tmp = A[n].split()        
                
        if len(tmp)==3:
            i = int(float(tmp[0])) - 1
            j = int(float(tmp[1])) - 1
            x = float(tmp[2])
        
            
            X.set(i,j, x*(1.0+0.0j))
            X.set(j,i, x*(1.0+0.0j))

            
    # Extract the sub-matrix of interest
    x_sub = CMATRIX(nstates1, nstates2)
    pop_submatrix(X, x_sub, act_sp1, act_sp2)
        
    return x_sub


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
    for i in range(0,nat):
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
    for i in range(0,nat):
        ln_indx = (nat+2)*md_iter1 + (i + 2)
        tmp = A[ln_indx].split()

        at = tmp[0]
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %s %15.8f  %15.8f  %15.8f\n" % (at, x, y, z)

    for i in range(0,nat):
        ln_indx = (nat+2)*md_iter2 + (i + 2)
        tmp = A[ln_indx].split()

        at = tmp[0]
        x = float(tmp[1])
        y = float(tmp[2])
        z = float(tmp[3])
        line = line + " %s %15.8f  %15.8f  %15.8f\n" % (at, x, y, z)

    return line 




def read_spectrum_unrestricted():
    """
    This function reads the eigenspectrum files
    produced by ErgoSCF in the case of spin-polarized
    (unrestricted) formulation.
    The names of the files are known, so the we do not need
    any arguments

    Note: the ordering of the orbitals in the ErgoSCF is:

    occupied_spectrum.txt       |    unoccupied_spectrum.txt
    ---------------------------------------------------------
    H                                       L
    H-1                                     L+1
    ...                                     ...
  
    So, the occupied orbitals order is reverted in the created array

    
    Args:  None

    Returns:
        (list, list): ( [e_a, e_b], [nocc, nvirt])
        
        * e_a ( list of floats ): energies of the alpha orbitals printed out
            this list contains both occupied and unoccupied orbitals
        * e_b ( list of floats ): energies of the beta orbitals printed out
            this list contains both occupied and unoccupied orbitals
        * nocc ( int ): the number of occupied orbitals, it should be less than 
            len(e_a) and len(e_b)
        * nvirt ( int ): the number of unoccupied orbitals, it should be less than 
            len(e_a) and len(e_b)

 
    """

    tmp = data_read.get_data_from_file2("occupied_spectrum_alpha.txt", [0])[0]
    nocc = len(tmp)
    e_a = []
    for i in range(nocc):
        e_a.append(tmp[nocc-1 - i])

    e_a = e_a + data_read.get_data_from_file2("unoccupied_spectrum_alpha.txt", [0])[0]
    nvirt = len(e_a) - nocc

    tmp = data_read.get_data_from_file2("occupied_spectrum_beta.txt", [0])[0]
    e_b = []
    for i in range(nocc):
        e_b.append(tmp[nocc-1 - i])

    e_b = e_b + data_read.get_data_from_file2("unoccupied_spectrum_beta.txt", [0])[0]


    return [e_a, e_b], [nocc, nvirt]



def read_mo_unrestricted(nocc, nvirt, act_space=None):
    """
    This function reads the MO files
    produced by ErgoSCF in the case of spin-polarized
    (unrestricted) formulation.
    The names of the files are known, so the we do not need
    to specify that 
    
    Args:  
        nocc ( int ): the number of occupied orbitals in the pool of 
           the printed out MO files
        nvirt ( int ): the number of unoccupied (virtual) orbitals in the pool of 
           the printed out MO files
        act_space ( list of ints ): indices of the orbitals we are interested in
            The indexing convention is relative to HOMO, that is 0 is HOMO, -1 is HOMO-1,
            1 is LUMO, etc.

    Returns:
        ( list ): ( [mo_a, mo_b] )
        
        * mo_a ( CMATRIX(nao, nocc+nvirt) ): alpha orbitals 
        * mo_b ( CMATRIX(nao, nocc+nvirt) ): beta orbitals
 
    """

    
    # Get dimensions:
    tmp_a = data_read.get_data_from_file2("homo_coefficient_vec_alpha.txt", [0])[0]
    nao = len(tmp_a)


    # Working active space
    ac = []     
    if act_space == None:
        ac = list(range(-(nocc-1), nvirt+1 ))
    else:
        ac = list(act_space)


    nmo =  len(ac)
    mo = [CMATRIX(nao, nmo), CMATRIX(nao, nmo)]


    for spin in [0, 1]:

        suff = ""
        if spin==0:
            suff = "alpha"
        elif spin==1:
            suff = "beta"

        # Mapping the relative (to HOMO) index of a MO to the
        # name of the file that contain this MO

        filenames = { 0 : F"homo_coefficient_vec_{suff}.txt",
                      1 : F"lumo_coefficient_vec_{suff}.txt"
                    }
        for j in range(1, nocc):
            filenames[-j] = F"occ_{j}_coefficient_vec_{suff}.txt"
        for j in range(2, nvirt+1):
            filenames[j] = F"unocc_{j-1}_coefficient_vec_{suff}.txt"

        # Now iterate over all entries in the active space
         
        for j in range(nmo):
            tmp_a = data_read.get_data_from_file2(filenames[ ac[j] ], [0])[0]
            for i in range(nao):
                mo[spin].set(i, j, tmp_a[i]*(1.0+0.0j))

    return mo





def read_spectrum_restricted():
    """
    This function reads the eigenspectrum files
    produced by ErgoSCF in the case of spin-unpolarized
    (restricted) formulation.
    The names of the files are known, so the we do not need
    any arguments

    Note: the ordering of the orbitals in the ErgoSCF is:

    occupied_spectrum.txt       |    unoccupied_spectrum.txt
    ---------------------------------------------------------
    H                                       L
    H-1                                     L+1
    ...                                     ...
  
    So, the occupied orbitals order is reverted in the created array
    
    Args:  None

    Returns:
        (list, list): ( [e_a], [nocc, nvirt])
        
        * e_a ( list of floats ): energies of the alpha orbitals printed out
            this list contains both occupied and unoccupied orbitals
        * nocc ( int ): the number of occupied orbitals, it should be less than 
            len(e_a) and len(e_b)
        * nvirt ( int ): the number of unoccupied orbitals, it should be less than 
            len(e_a) and len(e_b)

 
    """

    tmp = data_read.get_data_from_file2("occupied_spectrum.txt", [0])[0]
    nocc = len(tmp)

    e_a = []
    for i in range(nocc):
        e_a.append(tmp[nocc-1 - i])

    e_a = e_a + data_read.get_data_from_file2("unoccupied_spectrum.txt", [0])[0]
    nvirt = len(e_a) - nocc

    return [e_a], [nocc, nvirt]




def read_mo_restricted(nocc, nvirt, act_space=None):
    """
    This function reads the MO files
    produced by ErgoSCF in the case of spin-unpolarized
    (restricted) formulation.
    The names of the files are known, so the we do not need
    to specify that 
    
    Args:  
        nocc ( int ): the number of occupied orbitals in the pool of 
           the printed out MO files
        nvirt ( int ): the number of unoccupied (virtual) orbitals in the pool of 
           the printed out MO files
        act_space ( list of ints ): indices of the orbitals we are interested in
            The indexing convention is relative to HOMO, that is 0 is HOMO, -1 is HOMO-1,
            1 is LUMO, etc.


    Returns:
        ( list ): ( [mo_a] )
        
        * mo_a ( CMATRIX(nao, nocc+nvirt) ): alpha orbitals 
 
    """

    # Get dimensions:
    tmp_a = data_read.get_data_from_file2("homo_coefficient_vec.txt", [0])[0]
    nao = len(tmp_a)


    # Working active space
    ac = []     
    if act_space == None:
        ac = list(range(-(nocc-1), nvirt+1 ))
    else:
        ac = list(act_space)


    nmo =  len(ac)
    mo = [CMATRIX(nao, nmo)]


    for spin in [0]:

        suff = ""
        if spin==0:
            suff = "alpha"

        # Mapping the relative (to HOMO) index of a MO to the
        # name of the file that contain this MO

        filenames = { 0 : F"homo_coefficient_vec.txt",
                      1 : F"lumo_coefficient_vec.txt"
                    }
        for j in range(1, nocc):
            filenames[-j] = F"occ_{j}_coefficient_vec.txt"
        for j in range(2, nvirt+1):
            filenames[j] = F"unocc_{j-1}_coefficient_vec.txt"
        
        # Now iterate over all entries in the active space
         
        for j in range(nmo):
            tmp_a = data_read.get_data_from_file2(filenames[ ac[j] ], [0])[0]
            for i in range(nao):
                mo[spin].set(i, j, tmp_a[i]*(1.0+0.0j))

    return mo


def energies(en, nocc, nvirt, act_space=None):
    """
    This function creates a diagonal matrix of energies of the states
    included in the active space. 

    Args:  
        en ( list of floats ): energies of all (occupied + unoccupied) orbitals
            as read from the ErgoSCF
        nocc ( int ): the number of occupied orbitals in the pool of 
           the printed out MO files
        nvirt ( int ): the number of unoccupied (virtual) orbitals in the pool of 
           the printed out MO files
        act_space ( list of ints ): indices of the orbitals we are interested in
            The indexing convention is relative to HOMO, that is 0 is HOMO, -1 is HOMO-1,
            1 is LUMO, etc.

    Returns:
        ( CMATRIX(N,N) ): matrix with orbital energies (ordered from lowest to highest)
            Here, N is the size of the active space
        
    """

    # Working active space
    ac = []     
    if act_space == None:
        ac = list(range(-(nocc-1), nvirt+1 ))
    else:
        ac = list(act_space)

    nmo =  len(ac)
    E = CMATRIX(nmo, nmo)

    # Iterate over all entries in the active space
    for j in range(nmo):        
        E.set(j, j,  en[ ac[j] + (nocc-1) ] * (1.0+0.0j))

    return E
