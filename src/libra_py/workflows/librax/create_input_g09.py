#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#* Olga S. Bokareva
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \create_g09_input.py
#  This module defines the functions storing the template of G09 input file
#  and creating a G09 input file as JOB.inp.


import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def read_g09_inp_templ(inp_filename):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] inp_filename : The name of G09 input file used as template
    # This function returns 
    # templ : the format of G09 information before line with word "template".
    # 
    # Used in:  main.py/main

    f = open(inp_filename,"r")
    templ = f.readlines()
    f.close()

    N = len(templ)

    for i in range(0,N):
        s = templ[i].split()

        #if len(s) > 0 and s[0] == "template":
        if len(s) > 0 and s[0] == "Gfinput" and s[1]=="IOP(6/7=3)":
            ikeep = i # + 2
            break

    templ[ikeep+1:N] = []

    return templ
    

def write_g09_inp(label, Q, params, mol):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] label   A list of atomic labels (e.g. C,H,O)
    # \param[in] Q       A list of atomic charges (e.g. 6.0,1.0,8.0)
    # \param[in] params  A list of input parameters from g09_run.py 
    # \param[in] mol     An object containing nuclear DOF
    # This function returns the GAMESS input file (g09_inp)
    #
    # Used in:  md.py/run_MD

    g09_inp_templ = params["g09_inp_templ"]
    g09_inp = params["g09_inp"]
    mult = params["mult"]
    charge = params["charge"]

    g = open(g09_inp,"w")

    # Print the control section
    sz = len(g09_inp_templ)
    for i in xrange(sz):
        g.write(g09_inp_templ[i])

    # Important: Blank line
    g.write("\n")
    g.write("current calculation \n")
    g.write("\n")
    g.write("%d,%d \n" % (charge,mult))


    # Print coordinates
    Natoms = len(label)
    for i in xrange(Natoms):
        elt = label[i]
        q = Q[i]
        x = mol.q[3*i] 
        y = mol.q[3*i+1]
        z = mol.q[3*i+2]
        g.write("%s    %12.7f    %12.7f    %12.7f  \n"  % (elt, x, y, z) ) # q is not needed for g09

    g.write("\n\n")
    g.close()


    
