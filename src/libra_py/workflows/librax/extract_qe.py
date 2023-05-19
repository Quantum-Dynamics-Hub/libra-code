#*********************************************************************************
#* Copyright (C) 2016-2019 Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*^M
#*********************************************************************************/^M
"""
.. module:: extract_qe
   :platform: Unix, Windows
   :synopsis: This module implements functions for extracting data from the QE output file 
.. moduleauthor:: Ekadashi Pradhan, Alexey V. Akimov

"""

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *
import libra_py.packages.qe.methods as QE_methods
#>>>>>>>>>>>>>>>> UNCOMMENT THE SECTION BELOW, if THERE IS A PROBLEM WITH PATH
#cwd = "/projects/academic/alexeyak/ekadashi/libracode-dev/libracode-code/_build"
#print "Current working directory", cwd
#sys.path.insert(1,cwd+"/src/mmath")
#sys.path.insert(1,cwd+"/src/context")
#print "\nTest 1: Importing the library and its content"
#from libmmath import *
#from libcontext import *
#<<<<<<<<<<<<<<<<<<<<<<<<<



#def fermi_pop(e, nel, nspin, kT, el_st):
def excited_populations(e, nel, nspin, kT, el_st):
    """

    This function computes the occupation numbers based on Fermi populations with
    varied number of electrons, as described in our PCCP paper

    Args:
        e ( list of doubles ): Molecular orbital energies [ units: Ha]
        nel ( int ): The total number of electrons in the system
        nspin ( int ): The selector of the type of wavefunctions:

            * 1 - non-spin-polarized
            * 2 - spin-polarized calculations.

        kT  ( int ): the electronic energy (smearing constant) [ units: Ha]
        el_st ( int ): electronic state index: different fractional scheme for different electronic states.

    Returns:
        (list of doubles): occ_new: occupation numbers of the given MOs

    TODO:
        This implementation is crappy and it won't work in most cases, other than a few HOMO->LUMO transitions
        Needs to be generalized!

    """



    norbs = len(e)  # Total number of MOs in the active space    
    etol = 1e-10


    if el_st==0: # For S0 
        el_scheme = [0] # el_scheme is the fractional Fermi occupation scheme
                        # el_scheme is a list of integers representing indexes for orbitals
                        # with respect to the HOMO

    elif el_st==1: # For S1
        el_scheme = [-1,0,1]  # n_i^(N-1) + n_i^(N+1) - n_i^(N)

    elif el_st==2:  # For S2
        el_scheme = [-1,1,2]  # n_i^(N-1) + n_i^(N+2) - n_i^(N+1)

    else:
        print "Error: check out el_st, shouldn't be larger than 2 at this point.\n Exiting...\n"
        sys.exit(0)
    


    a = MATRIX(norbs,norbs)
    occ_tot = []

    for ia in el_scheme: # in el_scheme [-1,0,1], first element for N-1, second for N, and third for N+1 electrons.

        for i in xrange(norbs):
            a.set(i,i, e[i])

        if nspin == 2:                 # For spin-polarized calculations.
            Nel = nel/2 + nel%2 + ia   # Number of electrons in the alpha or beta spin orbital
            degen = 1.0                # One orbital can have 1 electrons, in this case of spin-polarization

        if nspin == 1:  # For non-polarized calculations
            Nel = nel + ia             # Total number of electrons.
            degen = 2.0                # One orbital can have 2 electrons, in case of non-polarized calculations


        bnds = order_bands(a)

        pop_fermi = populate_bands(Nel, degen, kT, etol, 1, bnds)  # 1 - use the fractional occ numbers


        tmp = []
        for item in pop_fermi:
            tmp.append(item[1])
        occ_tot.append(tmp) 

        #AVAL occ_tot.append([item[1] for item in pop_fermi]) 


    occ_new = []
    for ic in xrange(norbs):
        occ_new.append(0.0)

        if el_st==0: # For S0
            occ_new[ic] = occ_tot[0][ic] 

        elif el_st==1 or el_st==2: # For S1 and S2
            occ_new[ic] = occ_tot[0][ic] + occ_tot[2][ic] - occ_tot[1][ic]

        else:
            print "Warning: Higher schemes are not implemented yet, exiting"
            sys.exit(0)

    return occ_new




def gen_new_occ(ex_st,nel):
#def extract_orb_energy(HOMO):

    cwd1 = os.getcwd()
    en_alp = get_active_mo_en(cwd1+"/x%i.save/K00001/eigenval1.xml"%ex_st,nel)
    en_bet = get_active_mo_en(cwd1+"/x%i.save/K00001/eigenval2.xml"%ex_st,nel)
    #print "cwd",cwd1

    # push HOMO-1, HOMO, LUMO and LUMO+1 energies to get fermi energies
    # and fermi population, so, another function needed.
    occ_alp_new = fermi_pop(en_alp)

    # Similarly, compute fermi population for Beta spin, as this is a spin-polarized calculation
    occ_bet_new = fermi_pop(en_bet)

    return occ_alp_new, occ_bet_new


def write_qe_input(filename, nel,norb, flag_a, flag_b):
    #norb = 16
    nHOMO = nel/2 + nel%2 # HOMO orbital index
    if norb%10 ==0:
        nl_spin_orb = norb/10      # Number of lines in spin orbital
    else:
        nl_spin_orb = norb/10 + 1  # For Ethylene, it is 2, 12/10 + 1

    #print "nHOMO=",nHOMO
    fa = filename.split('.')
    fb = fa[0]+".scf_wrk.in"
    f = open(fb,"r+")
    a = f.readlines()
    N = len(a)
    f.close()
    f = open(fb, "w")
    i_homo = -1 # initializing with random value
    i_homo_bet = -1 # initializing
    homo_idx = (nHOMO%10) -1
    for i in range(0,N):
        s = a[i].split()
        if len(s)>0 and s[0] =="OCCUPATIONS":
            i_alp = i+1
            #print "i_alp=",i_alp
            i_homo = i_alp + nHOMO/10
            i_bet = i + nl_spin_orb + 2
            i_homo_bet = i_bet + nHOMO/10
            #print "i_homo=",i_homo

        if i==i_homo:
            s[4],s[5],s[6],s[7] = flag_a[0][1],flag_a[1][1],flag_a[2][1],flag_a[3][1]
            a[i] = " ".join(str(x) for x in s)+'\n'
        if i==i_homo_bet:
            s[4],s[5],s[6],s[7] = flag_b[0][1],flag_b[1][1],flag_b[2][1],flag_b[3][1]
            a[i] = " ".join(str(x) for x in s)+'\n'

        f.write(a[i])
    f.close()


def robust_cal_extract(filename, ex_st, nel, flag_a, flag_b):
    write_qe_input(filename, nel,flag_a, flag_b)
    exe_espresso(ex_st)
    tot_ene,iforce = extract_ene_force(filename)
    return tot_ene, iforce



def qe_extract(ex_st, active_space, nspin):
    """
 
    This function reads the extracting Quantum Espresso files produced in
    x0.save, x1.save, ...  and x0.export, x1.export, ... folders

    Args:
        ex_st ( int ): index of excited state to read info about, this index should 
            be consistent with the naming of the output directories, which should be 
            labeled according to prefixes "x0", "x1", etc.
        active_space ( list of ints ): Indices (starting from 1) of the MOs to include in
            calculations (and to read from the QE output files) 
        nspin ( int ): type of spin-polarization in the calculations:
  
            * 1: non-spin-plarized
            * 2: spin-plarized

    Returns:
        tuple: (info, MO_a, MO_b), where:

            * info ( dictionary ): seealso: ```QE_methods.read_qe_schema``` for mor details
            * MO_a ( CMATRIX(npw, norbs) ): alpha-spin orbitals, where npw is the number of
                plane-waves in the basis, norbs - is the number of orbitals included in the
                active space norbs = len(active_space)
            * MO_b ( CMATRIX(npw, norbs) ): beta-spin orbitals, same as ```MO_a```

    """

    info = QE_methods.read_qe_schema("x%i.save/data-file-schema.xml" % (ex_st), 0)

    MO_a, MO_b = None, None
    if nspin <= 1:
        # Read the wavefunctions:
        MO_a = QE_methods.read_qe_wfc("x%i.export/wfc.1" % (ex_st), active_space, 0)
        MO_b = CMATRIX(MO_a)

    if nspin == 2:
        # Read the wavefunctions:
        MO_a = QE_methods.read_qe_wfc("x%i.export/wfc.1" % (ex_st), active_space, 0)
        MO_b = QE_methods.read_qe_wfc("x%i.export/wfc.2" % (ex_st), active_space, 0)

    return info, MO_a, MO_b

