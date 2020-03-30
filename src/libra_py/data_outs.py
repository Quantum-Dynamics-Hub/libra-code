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
.. module:: data_outs
   :platform: Unix, Windows
   :synopsis: 
       This module implements various functions for printing out data

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


def print_matrix(x):

    if type(x).__name__ in ["MATRIX", "IMATRIX", "CMATRIX"]:

        nrows = x.num_of_rows
        ncols = x.num_of_cols

        for row in range(nrows):
            line = ""
            for col in range(ncols):
                line = line + F"{x.get(row, col)}  "
            print(line)

    elif type(x).__name__ == "MATRIX3x3":
        print(F" xx = {x.xx}  xy = {x.xy}  xz = {x.xz}")
        print(F" yx = {x.yx}  yy = {x.yy}  yz = {x.yz}")
        print(F" zx = {x.zx}  zy = {x.zy}  zz = {x.zz}")



def show_matrix_pyplot(X, set_diag_to_zero=0):
    """

    This function prints a matrix in a file in the format
    recognized by the matplotlib for 2D surfaces and maps

    Args:
        X ( MATRIX(N,M) ): object containing the data to be printed out
        set_diag_to_zero ( int ): wheather print out the diagonal elements as 0.0 

            * 0 - no, print out the diagonal elements as they are [ default ]
            * 1 - yes, print out the diagonal elements as 0.0

    Returns: 
        tuple: (x, y, z):
        
            * x ( list of doubles ): x grid 
            * y ( list of doubles ): y grid 
            * z ( list of doubles ): values of z = z(x,y)

    """

    ncol, nrow = X.num_of_cols, X.num_of_rows
    x, y = [], []
    for i in range(0,ncol):
        x.append(i)
    for i in range(0,nrow):
        y.append(i)

    
    z = []    
    for i in range(0,nrow):
        z_i = []
        for j in range(0,ncol):
            val = X.get(i,j)
            if set_diag_to_zero==1:
                if i==j:
                    val = 0.0
            z_i.append(val)
        z.append(z_i)
         
    return x, y, z



def show_matrix_splot(X, filename, set_diag_to_zero=0):
    """

    This function prints a matrix in a file in the format
    recognized by the gnuplot's splot function (for 2D surfaces and maps)

    Args:
        X ( MATRIX(N,M) ): object containing the data to be printed out
        filename ( string ): the name of the file to where the data will be printed out
        set_diag_to_zero ( int ): wheathere print out the diagonal elements as 0.0 

            * 0 - no, print out the diagonal elements as they are [ default ]
            * 1 - no, print out the diagonal elements as 0.0

    Returns:
        string: the line containing the information also printed out to the file
            This string is in a format suitable for 2D plotting with gnuplot

    """

    ncol, nrow = X.num_of_cols, X.num_of_rows

    line = ""
    for i in range(0,nrow):
        for j in range(0,ncol):
            val = X.get(i,j)
            if set_diag_to_zero==1:
                if i==j:
                    val = 0.0
            line = line + "%4i %4i %8.5f \n" % (i, j, val)
        line = line + "\n"
        
    f = open(filename, "w")
    f.write(line)
    f.close()
 
    return line



def printout(t, pops, Hvib, outfile):
    """
    t - time [a.u.] 
    pops - [MATRIX] - populations
    Hvib - [CMATRIX] - vibronic Hamiltonian
    outfile - filename where we'll print everything out
    """

    nstates = Hvib.num_of_cols


    line = "%8.5f " % (t)
    P, E = 0.0, 0.0
    for state in range(0,nstates):
        p, e = pops.get(state,0), Hvib.get(state, state).real
        P += p
        E += p*e
        line = line + " %8.5f  %8.5f " % (p, e)
    line = line + " %8.5f  %8.5f \n" % (P, E)

    f = open(outfile, "a") 
    f.write(line)
    f.close()


    

def add_printout(i, pop, filename):
    # pop - CMATRIX(nstates, 1)

    f = open(filename,"a")
    line = "step= %4i " % i    

    tot_pop = 0.0
    for st in range(0,pop.num_of_cols):
        pop_o = pop.get(st,st).real
        tot_pop = tot_pop + pop_o
        line = line + " P(%4i)= %8.5f " % (st, pop_o)
    line = line + " Total= %8.5f \n" % (tot_pop)
    f.write(line)
    f.close()



def list2bin(x, name):
    """
    Converts a list of doubles into binary file    
    
    """
    sz = len(x)
    X = MATRIX(sz,1)
    
    for i in range(sz):
        X.set(i,0, x[i])
    
    X.bin_dump(name+"_%i" % (sz))


def bin2list(name, sz):
    """
    Read in a binary file into a matrix and then
    convert it to a Python list
    
    """
    
    X = MATRIX(sz, 1)
    X.bin_load(name+"_%i" % (sz))
        
    x = []
    for i in range(sz):
        x.append(X.get(i,0))
    
    return x

