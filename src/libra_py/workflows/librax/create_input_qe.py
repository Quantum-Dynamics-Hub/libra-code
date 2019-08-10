#*********************************************************************************
#* Copyright (C) 2016-2019 Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: create_input_qe
   :platform: Unix, Windows
   :synopsis: This module implements the functions which prepare QE input with different occupation schemes
.. moduleauthor:: Ekadashi Pradhan, Alexey V. Akimov

"""

import os
import sys
import math
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

# TODO: Get rid of import *
#from libra_py import *
import util.libutil as comn


def read_qe_inp_templ(inp_filename):
    """

    Reading and storing template input file for QE calculations. The input file is essentially a
    normal input file, but we store only the constant section (control option), not the
    coordinates. The latter will be updated at each iteration using the propagated objects

    Args:
        inp_filename ( string ): The name of the initial input file, which will serve as a template


    """

    f = open(inp_filename,"r")
    templ = f.readlines()
    f.close()

    for a in templ:
        s = a.split()
        if len(s) > 0 and s[0] == "celldm(1)" and s[1] == "=":
            sa = s[2].split(',')
            cell_dm = float(sa[0])
            break

    # Find the line preceeding the actual atomic coordinates
    for a in templ:
        s = a.split()
        if len(s) > 0 and s[0] == "ATOMIC_POSITIONS":
            ikeep = templ.index(a)
            break

    N = len(templ)
    # Blank space for the atomic positions
    templ[ikeep+1:N] = []
    for i in xrange(ikeep+1):
        print templ[i]


    return  templ



def excitation_to_qe_occ(norb, nel, state):
    """

    This function converts Libra "excitation" objects to the QE occupation scheme

    Args:
        norb ( int ): the number of KS orbitals to consider
        nel ( int ): the number of electrons to distribute among these orbitals
        state ( Libra.excitation object ): the excitation representing Slater determinant

    Returns:
        (list, list, list): (occ, occ_alp, occ_bet):

            * occ ( list of doubles ): total occupation of the orbitals (2 is max, for non-spin-polarized)
            * occ_alp ( list of doubles ): occupations of alpha orbitals (for spin-polarized)
            * occ_bet ( list of doubles ): occupations of beta orbitals (for spin-polarized)

    """

    # Number of occupied alpha and beta orbitals
    nocc_alp = nel/2  # integer division!
    nocc_bet = nel - nocc_alp
    homo = nocc_alp -1  # changed indices in order to accomodate with python numbering

    # Generate reference (ground state) occupation scheme for alpha and beta orbitals
    gs_alp = []
    gs_bet = []
    for i in xrange(norb):
        if i<nocc_alp:
            gs_alp.append([i,1.0])
        else:
            gs_alp.append([i,0.0])

        if i<nocc_bet:
            gs_bet.append([i,1.0])
        else:
            gs_bet.append([i,0.0])

    # Compute indices of the orbitals involved in the excitation
    a = state.from_orbit[0] + homo  
    a_s = state.from_spin[0]  # +1 for alp, -1 for bet
    b = state.to_orbit[0] + homo   
    b_s = state.to_spin[0]    # +1 for alp, -1 for bet

    # Do separate alpha and beta excitations
    # Here we use "excite" function from the core of Libra package
    ex_alp = []
    ex_bet = []
    if a_s==1 and a_s == b_s:
        ex_alp = excite(a,b,gs_alp) 
        ex_bet = gs_bet
    elif a_s==-1 and a_s == b_s:
        ex_alp = gs_alp
        ex_bet = excite(a,b,gs_bet)
    else:
        print "Error in qe_occ: excitations with spin flip are not allowed yet\nExiting...\n"
        sys.exit(0)


    occ = [0.0]*norb
    occ_alp = [0.0]*norb
    occ_bet = [0.0]*norb

    for i in xrange(norb):
        occ_alp[i] = ex_alp[i][1]
        occ_bet[i] = ex_bet[i][1]
        occ[i] = ex_alp[i][1] + ex_bet[i][1]

    return occ, occ_alp, occ_bet



def print_occupations(occ):
    """

    This function transforms the list of occupations into a formatted
    text. The format is consistent with the QE input

    Args:
        occ ( list of floats ): Occupation numbers 

    Returns:
        string: line: a multiline text that contains the occupation numbers
            suitable for QE input
    """
    
    #line = "OCCUPATIONS\n"
    line = ""
    count = 0
    for f in occ:
        line = line + "%15.12f " % f 
        count = count +1
        if count % 10 ==0:
            line = line + "\n"
    line = line + "\n"

    return line



def write_qe_input(ex_st, label, mol, params,occ,occ_alp,occ_bet,restart_flag):
##
# Creates the Quantum Espresso input using the data provided
# \param[in] ex_st The index of the excited state we want to compute - so it controls which input file
# to create and how to do this
# \param[in] label Element symbols for all atoms in the system (list of strings)
# \param[in] mol The object containing nuclear DOF
# \param[in] params The general control parameters (dictionary)
# \param[in] occ a list of occupation numbers of the total orbitals (doubly degenerate)
# \param[in] occ_alp a list of occupation numbers of the alpha orbitals
# \param[in] occ_bet a list of occupation numbers of the beta orbitals
# \param[in] restart_flag index 10 is used at this point for the very first iteration, 11 for consecutive
#            steps where restart from the previous wavefunction and density is performed.

    HOMO = params["nel"]/2 - 1 # It must be integer, This is HOMO index
    excitation = params["excitations"][ex_st]
    qe_inp = "x%i.scf_wrk.in" % ex_st
    nspin = params["nspin"]

    qe_inp_templ = params["qe_inp_templ"][ex_st]
    cell_dm = params["alat"]
    pp = qe_inp.split('.')
    prefix = pp[0]
    f = open(qe_inp, "w")     

    # Write control parameters section
    for a in qe_inp_templ:
        aa = a.split()

        # Make sure the content of input files is consistent with the file names
        if len(aa) > 0:
            if aa[0] == "prefix":
                a = "  prefix = '%s',\n" % (prefix)

        # Force the calculations to re-use existing files
        if len(aa) > 0:
            if aa[0] == "&ELECTRONS" and restart_flag==11:
                a = "&ELECTRONS \n"
                a = a + " startingwfc = 'file', \n"
                a = a + " startingpot = 'file', \n"

        if len(aa) > 0:
            if aa[0] == "&ELECTRONS" and restart_flag==10:
                a = "&ELECTRONS \n"

        # Change the number of max SCF iterations in the input file
        if len(aa) > 0:
            if aa[0] == "electron_maxstep" and restart_flag > 9:
                a = " electron_maxstep = %i, \n " % (params["scf_itr"])


        f.write(a)
    f.write("\n")

    # Write atom name and coordinatess
    Natoms = len(label)
    B_to_A = 1.0/cell_dm   # Bohr to Angstrom conversion

    for k in xrange(Natoms):
        atms = label[k]
        x = B_to_A*mol.q[3*k]
        y = B_to_A*mol.q[3*k+1]
        z = B_to_A*mol.q[3*k+2]
        f.write("%s    %12.7f    %12.7f    %12.7f  \n"  % (atms, x, y, z) )


    ####################################################
    #  Give conditional statement, if 
    #  if flag1 == -1:
    #      generate occupation number using fermi population
    #  elif flag1 == 0:
    #      continue the general sequence of writing input file       
    ####################################################

    # Write occupation
    f.write("\nOCCUPATIONS\n")

    # Single excitations with no spin-polarization
    if nspin <= 1:
        f.write(print_occupations(occ))
        f.write("\n")

    # Single excitations with spin-polarization 
    if nspin >1:
        f.write(print_occupations(occ_alp))
        f.write("\n")
        f.write(print_occupations(occ_bet))        
    f.close()



def write_qe_input_first(filename1, filename2, occ, occ_alp, occ_bet, nspin, scf_iter, restart_flag):
    """

    Creates the QE input file for the very first step

    Args:
    filename1 ( string ): the QE input filename that is used as a reference (template)
    filename2 ( string ): the QE input filename that will be produced
    occ ( list of doubles ): occupation numbers of the orbitals in non-polarized case (doubly degenerate)
    occ_alp ( list of doubles ): occupation numbers of the alpha orbitals (in spin-polarized case)
    occ_bet ( list of doubles ): occupation numbers of the beta orbitals (in spin-polarized case)
    nspin ( int ): selection of the spin-polarization method:

        * 1: spin restricted (non-polarized)
        * 2: spin unrestricted (polarized)

    scf_iter ( int ): The number of SCF iterations to do in the QE calculations
    restart_flag ( int ): 
        * 0: don't change the number of max SCF iterations
        * 1: DO change the number of max SCF iterations, but don't use the previous wfc and potential
        * 2: DO change the number of max SCF iterations, and DO use the previous wfc and potential

    """
 
    f = open(filename2,"r")
    a = f.readlines()
    f.close()

    N = len(a)  # number of lines

    # Now re-write the file
    i_alp = 0
    for i in range(0,N):
        s = a[i].split()
        if len(s)>0 and s[0] =="OCCUPATIONS":
            i_alp = i+1
    a[i_alp:N] = []


    f = open(filename2, "w")
    for i in range(0,i_alp):
        aa = a[i].split()
        if len(aa) >0 and aa[0] == "&ELECTRONS" and restart_flag==2:
            a[i] = "&ELECTRONS \n startingwfc = 'file', \n startingpot = 'file', \n"
        if len(aa) >0 and aa[0] == "&ELECTRONS" and restart_flag==1:
            a[i] = "&ELECTRONS \n"
        if len(aa) >0 and aa[0] == "electron_maxstep" and restart_flag > 0:
            #a = " electron_maxstep = 2, \n "
            a[i] = " electron_maxstep = %i, \n " % (scf_iter) 
        f.write(a[i])

    # Write occupation
    # Single excitations with no spin-polarization
    if nspin <= 1:
        f.write(print_occupations(occ))
        f.write("\n")

    # Single excitations with spin-polarization 
    if nspin > 1:
        f.write(print_occupations(occ_alp))
        f.write("\n")
        f.write(print_occupations(occ_bet))

    f.close()














