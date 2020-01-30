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
.. module:: dynamics_io
   :platform: Unix, Windows
   :synopsis: This module implements the read/write functions specifically designed to work with the dynamics module
       List of functions:

           * add_intlist2file(filename, t, X)
           * add_doublelist2file(filename, t, X)
           * add_matrix2file(filename, t, X)
           * add_cmatrix2file(filename, t, X)
           * file2intlist(filename)
           * file2doublelist(filename)
           * file2matrix(filename, nrows, ncols)
           * file2matrix(filename, nrows, ncols)
           * print_results12(i, dt, res, prefix, file_output_level)
           * print_results3(i, dt, res, prefix, file_output_level, tr)
           * read_results(prefix, file_output_level, nadi, ndia, ndof, ntraj)

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




def print_results12(i, dt, res, prefix, file_output_level):
    """
    Print out the properties at given timestep i to the files

    Args:
        i ( int ): index of the timestep
        dt ( double ): time step size [ units: a.u. ]
        res ( tuple ): contains the properties (at given time) to print out, the order is as follows:

            * q = res[0] ( MATRIX(ndof, ntraj) ): coordiantes of all trajectories
            * p = res[1] ( MATRIX(ndof, ntraj) ): momenta of all trajectories
            * Ekin = res[2] ( double ): average kinetic energy of the ensemble
            * Epot = res[3] ( double ): average potential energy of the ensemble
            * Etot = res[4] ( double ): average total energy of the ensemble
            * dEkin = res[5] ( double ): fluctuation of the average kinetic energy of the ensemble
            * dEpot = res[6] ( double ): fluctuation of the average potential energy of the ensemble
            * dEtot = res[7] ( double ): fluctuation of the average total energy of the ensemble
            * Cadi = res[8] ( CMATRIX(nadi, ntraj) ): amplitudes of adiabatic states for all trajectories
            * Cdia = res[9] ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states for all trajectories
            * dm_adi = res[10] ( CMATRIX(nadi, nadi) ): adiabatic density matrix averaged over all trajectories
            * dm_dia = res[11] ( CMATRIX(ndia, ndia) ): diabatic density matrix averaged over all trajectories
            * pops = res[12] ( list of nadi doubles ): average populations of all adiabatic states
            * states = res[13] ( list of ntraj ints ): electronic states of all trajectories

        prefix ( string ): the common prefix of the output (usually the directory name)
        file_output_level ( int ): the verbosity of the output

    Returns:
        None

    """

    q, p = res[0], res[1]
    Ekin, Epot, Etot = res[2], res[3], res[4]
    dEkin, dEpot, dEtot = res[5], res[6], res[7]
    Cadi, Cdia = res[8], res[9]
    dm_adi, dm_dia = res[10], res[11]
    pops = res[12]
    states = res[13]
    
    # File output
    if file_output_level >= 1:
        add_doublelist2file("%s/energies.txt" % (prefix), i*dt, [Ekin, Epot, Etot, dEkin, dEpot, dEtot] )
        add_cmatrix2file("%s/D_adi.txt" % (prefix), i*dt, dm_adi )
        add_cmatrix2file("%s/D_dia.txt" % (prefix), i*dt, dm_dia )
        add_matrix2file("%s/SH_pop.txt" % (prefix), i*dt, pops )

    # File output
    if file_output_level >= 2:        
        add_matrix2file("%s/q.txt" % (prefix), i*dt, q )
        add_matrix2file("%s/p.txt" % (prefix), i*dt, p )
        add_cmatrix2file("%s/C_adi.txt" % (prefix), i*dt, Cadi )
        add_cmatrix2file("%s/C_dia.txt" % (prefix), i*dt, Cdia )
        add_intlist2file("%s/states.txt" % (prefix), i*dt, states )        



def print_results3(i, dt, res, prefix, file_output_level, tr):
    """
    Print out the properties at given timestep i to the files

    Args:
        i ( int ): index of the timestep
        dt ( double ): time step size [ units: a.u. ]
        res ( tuple ): contains the properties (at given time) to print out, the order is as follows:

            * hvib_adi = res[0] ( CMATRIX(nadi, nadi) ): adiabatic Hamiltonian for a given trajectory
            * hvib_dia = res[1] ( CMATRIX(ndia, ndia) ): diabatic Hamiltonian for a given trajectory
            * St = res[2] ( CMATRIX(nadi, nadi) ): time-overlap of the adiabatic states for a given trajectory
            * U = res[3] ( CMATRIX(ndia, nadi) ): dia-to-adi transformation matrix for a given trajectory

        prefix ( string ): the common prefix of the output (usually the directory name)
        file_output_level ( int ): the verbosity of the output
        tr ( int ): the index of the trajectory that we handle

    Returns:
        None

    """

    hvib_adi = res[0]
    hvib_dia = res[1]
    St = res[2]
    U = res[3]

    # File output
    if file_output_level >= 3:
        add_cmatrix2file("%s/Hvib_adi_%i.txt" % (prefix, tr), i*dt, hvib_adi )
        add_cmatrix2file("%s/Hvib_dia_%i.txt" % (prefix, tr), i*dt, hvib_dia )
        add_cmatrix2file("%s/St_%i.txt" % (prefix, tr), i*dt, St)
        add_cmatrix2file("%s/basis_transform_%i.txt" % (prefix, tr), i*dt, U )



def read_results(prefix, file_output_level, nadi, ndia, ndof, ntraj):
    """
    This function 
    """
    
    obs_T = None
    obs_q, obs_p = None, None
    obs_Ekin, obs_Epot, obs_Etot = None, None, None
    obs_dEkin, obs_dEpot, obs_dEtot = None, None, None
    obs_Cadi, obs_Cdia = None, None 
    obs_dm_adi, obs_dm_dia = None, None
    obs_pop, obs_states = None, None
    obs_hvib_adi, obs_hvib_dia, obs_St, obs_U = None, None, None, None

    
    if file_output_level >= 1:
        res = data_read.get_data_from_file2("%s/energies.txt" % (prefix), range(0, 7) )
        obs_T = res[0]
        obs_Ekin = res[1]
        obs_Epot = res[2]
        obs_Etot = res[3]
        obs_dEkin = res[4]
        obs_dEpot = res[5]
        obs_dEtot = res[6]
                
        obs_dm_adi = file2cmatrix("%s/D_adi.txt" % (prefix), nadi, nadi)
        obs_dm_dia = file2cmatrix("%s/D_dia.txt" % (prefix), ndia, ndia)
        obs_pop = file2matrix("%s/SH_pop.txt" % (prefix), nadi, 1)
        
    if file_output_level >= 2:
        obs_q = file2matrix("%s/q.txt" % (prefix), ndof, ntraj)
        obs_p = file2matrix("%s/p.txt" % (prefix), ndof, ntraj)
        obs_C_adi = file2cmatrix("%s/C_adi.txt" % (prefix), nadi, ntraj)
        obs_C_dia = file2cmatrix("%s/C_dia.txt" % (prefix), ndia, ntraj)
        obs_states = file2intlist("%s/states.txt" % (prefix))       
        
    if file_output_level >= 3:
        obs_hvib_adi, obs_hvib_dia, obs_St, obs_U = [], [], [], []
        for tr in range(ntraj):                        
            hvib_adi = file2cmatrix("%s/Hvib_adi_%i.txt" % (prefix, tr), nadi, nadi)
            hvib_dia = file2cmatrix("%s/Hvib_dia_%i.txt" % (prefix, tr), nadi, nadi)
            St = file2cmatrix("%s/St_%i.txt" % (prefix, tr), nadi, nadi)
            U = file2cmatrix("%s/basis_transform_%i.txt" % (prefix, tr), ndia, nadi)
            
            obs_hvib_adi.append(hvib_adi)
            obs_hvib_dia.append(hvib_dia)
            obs_St.append(St)
            obs_U.append(U)                


    return obs_T, obs_q, obs_p, obs_Ekin, obs_Epot, obs_Etot, obs_dEkin, obs_dEpot, obs_dEtot, \
           obs_Cadi, obs_Cdia, obs_dm_adi, obs_dm_dia, obs_pop, obs_states, obs_hvib_adi, obs_hvib_dia, obs_St, obs_U

