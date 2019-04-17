#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Kosuke Sato, Alexey V. Akimov
#* 
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file gamess_to_libra.py 
# This module implements the functions that extract parameters from the gamess output file:
# atomic forces , molecular energies, molecular orbitals, and atomic basis information.
# The forces are used for simulating Classical MD on Libra 
# and the others for calculating time-averaged energies and Non-Adiabatic Couplings(NACs).

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from extract_qe import *
from overlap import *
from hamiltonian_el import *
from create_input_qe import *
from misc import *
from spin_indx import *



def exe_espresso(i):
##
# Function for executing calculations using Quantum Espresso
# once the calculations are finished, all the temporary data are
# deleted
# \param[in] inp The name of the input file
# \param[in] out The name of the output file
#
    inp = "x%i.scf_wrk.in" % i # e.g. "x0.scf_wrk.in"
    out = "x%i.scf.out" % i    # e.g. "x0.scf.out"
    inexp = "x%i.exp.in" % i   # e.g. "x0.exp.in"
    outexp = "x%i.exp.out" % i # e.g "x0.exp.out"

    os.system("srun pw.x < %s > %s" % (inp,out))
    os.system("srun pw_export.x < %s > %s" % (inexp,outexp))

    # Delete scratch directory and unecessary files
    #os.system("rm *.dat *.wfc* *.igk* *.mix*")
    #os.system("rm -r *.save") # not sure if we  need to remove this directory



def qe_to_libra(params, E, sd_basis, label, mol, suff, active_space):
    ## 
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params :  contains input parameters , in the directory form
    # \param[in, out] E  :  molecular energies at "t" old, will be updated (MATRIX object)
    # \param[in] sd_basis :  basis of Slater determinants at "t" old (list of CMATRIX object)
    # \param[in] label : the list of atomic names 
    # \param[in] mol : the object of Nuclear type - contains the info about molecular geometry
    # \param[in] suff : The suffix to add to the name of the output files
    # this suffix is now considered to be of a string type - so you can actually encode both the
    # iteration number (MD timestep), the nuclear cofiguration (e.g. trajectory), and any other
    # related information
    # \param[in] active_space The list of indices (starting from 1) of the MOs to include in
    # calculations (and to read from the QE output files)

    #
    # This function outputs the files for excited electron dynamics
    # in "res" directory.
    # It returns the forces which act on the atoms.
    # Also, it returns new atomic orbitals, molecular energies, and
    # molecular coefficients used for calculating time-averaged
    # molecular energies and Non-Adiabatic Couplings(NACs).
    #
    # Used in: md.py/run_MD

    nstates = len(params["excitations"])

    sd_basis2 = SDList()    # this is a list of SD objects. Eeach represents a Slater Determinant
    all_grads = [] # this will be a list of lists of VECTOR objects
    E2 = MATRIX(nstates,nstates)
    nspin = params["nspin"]
    nel = params["nel"]
    HOMO = nel/2 + nel%2 -1
    
    #======== Run QE calculations and get the info at time step t+dt ========
    
    for ex_st in xrange(nstates): # for each excited configuration
                                  # run a separate set of QE calculations
        #idx = params["excitations"][ex_st]
        ###########################################
        #write_qe_input(ex_st,label,mol,params)
        #exe_espresso(ex_st)
        #nspin = params["nspin"]
        #tot_ene, label, R, grads, mo_pool_alp, mo_pool_bet, norb, nel, nat, alat = qe_extract("x%i.scf.out" % ex_st, flag, active_space, ex_st,nspin)
        ###########################################
        flag = 0
        mx_itr = params["max_iteration"]
        #while flag1! = 0: #for i in xrange(5):
        #    write_qe_input(ex_st,label,mol,params,flag1,occ_a,occ_b)
        #    exe_espresso(ex_st)

            #flag = 0
        #    nspin = params["nspin"]
            # Return flag1, if this is -1, then continue the loop, it != -1, then break the loop
        #    flag1, tot_ene, label, R, grads, mo_pool_alp, mo_pool_bet, norb, nel, nat, alat, occ_a, occ_b = qe_extract("x%i.scf.out" % ex_st, flag, active_space, ex_st,nspin)

        #occ_alp,occ_bet = excitation_occ(params)
        excitation = params["excitations"][ex_st]
        occ, occ_alp, occ_bet = excitation_to_qe_occ(params, excitation)
        status = -1
        restart_flag = 0
        coount = 0
        while status != 0: #for i in xrange(5):
            coount = coount + 1
            write_qe_input(ex_st,label,mol,params,occ,occ_alp,occ_bet,restart_flag)
            exe_espresso(ex_st)
            status = check_convergence("x%i.scf.out" % ex_st) # returns 0 if SCF converges, 1 if not converges
            if status == 0:
                tot_ene, label, R, grads, mo_pool_alp, mo_pool_bet, norb, nel, nat, alat = qe_extract("x%i.scf.out" % ex_st, active_space, ex_st, nspin, flag)

            else:
                if coount==1:
                    restart_flag = 10
                if coount>mx_itr:  # Maximum 30 iteration, else exit
                    print "Warning! Maximum iteration for the Fermi scheme reached, please check electronic_smearing parameter - exiting"
                    sys.exit(0)
                else:
                    restart_flag = 11
                if params["nspin"] == 2:
                    en_alp = qe_extract_eigenvalues("x%i.save/K00001/eigenval1.xml"%ex_st,nel)
                    en_bet = qe_extract_eigenvalues("x%i.save/K00001/eigenval2.xml"%ex_st,nel)
                    occ_alp = fermi_pop(en_alp,nel,params["nspin"],params["electronic_smearing"],ex_st)
                    occ_bet = fermi_pop(en_bet,nel,params["nspin"],params["electronic_smearing"],0)
                if params["nspin"] == 1:
                    en_orb = qe_extract_eigenvalues("x%i.save/K00001/eigenval.xml"%ex_st,nel)
                    occ = fermi_pop(en_orb,nel,params["nspin"],params["electronic_smearing"])

#            if flag1!=-1: # when iforce != -1, then break this loop and continue
#                break

        #for la in xrange(5):
        #    a1,b1 = gen_new_occ(ex_st) # This can be inserted inside qe_extract such that is gives a1 and b1 as output. i.e., extracts additional information
        #    Then it is easy to just write input and rerun calculation
        #    tot_ene, iforce = robust_cal_extract(filename, ex_st, nel, a1,b1)
        #    if iforce !=-1:
        #        break


        #mo_pool_alp = CMATRIX(mo_pool) 
        #mo_pool_bet = CMATRIX(mo_pool)


        homo = nel/2 +  nel % 2
        alp,bet = index_spin(params["excitations"][ex_st], active_space, homo) 

        # use excitation object to create proper SD object for different excited state
        sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet) )
        sd_basis2.append(sd)
        all_grads.append(grads)
        
        E2.set(ex_st, ex_st, tot_ene)

    # calculate overlap matrix of Slater determinant basis states
    P11 = SD_overlap(sd_basis,  sd_basis)
    P22 = SD_overlap(sd_basis2, sd_basis2)
    P12 = SD_overlap(sd_basis,  sd_basis2)
    P21 = SD_overlap(sd_basis2, sd_basis)




    ### TO DO: In the following section, we need to avoid computing NAC matrices in the full
    # basis. We will need the information on cropping, in order to avoid computations that 
    # we do not need (the results are discarded anyways)
    # calculate molecular energies and Non-Adiabatic Couplings(NACs) on MO basis

    # ********************************************************************************************
    # The names "E_mol" and "D_mol" are confusing when you have already constructed Energy and NACs           
    # based on SD basis sets.
    # ********************************************************************************************
    #S_mol = average_S(P11,P22)  # Average S(t+dt/2) = (S(t) + S(t+dt))/2.0
    S_mol = 0.50 * (P11 + P22)  # Seperate module is not required, directly averaged here
    E_mol = average_E(E,E2)

    if params["non-orth"] ==1:
        D_mol = P12
    else:
        D_mol = NAC(P12,P21,params["dt_nucl"])

    # reduce the matrix size
    #E_mol_red = reduce_matrix(E_mol,params["min_shift"], params["max_shift"],params["HOMO"])
    #D_mol_red = reduce_matrix(D_mol,params["min_shift"], params["max_shift"],params["HOMO"])
    ### END TO DO

    if params["print_mo_ham"]==1:
        E_mol.show_matrix(params["mo_ham"]+"full_re_Ham_"+suff)
        D_mol.show_matrix(params["mo_ham"]+"full_im_Ham_"+suff)
        #E_mol_red.show_matrix(params["mo_ham"] + "reduced_re_Ham_" + suff)
        #D_mol_red.show_matrix(params["mo_ham"] + "reduced_im_Ham_" + suff)



    # store "t+dt"(new) parameters on "t"(old) ones
    E = MATRIX(E2)  # update energy
                    # the returned energy E_mol is at t+dt/2

    # Returned data:
    # 
    # Grad: Grad[k][i] - i-th projection of the gradient w.r.t. to k-th nucleus (i = 0, 1, 2)
    # data: a dictionary containing transition dipole moments
    # E_mol: the matrix of N-el orbital (total) energies at t+dt/2 in the reduced (active) space
    # D_mol: the matrix of the NACs computed with SD orbitals at t+dt/2 in the reduced (active) space
    # E2: the matrix of the N-el energies at t+dt (present state)
    # sd_basis2: the list of reduced SD (active space orbitals), representing all computed states
    # all_grads: the gradients on all atoms for all excited states, such that all_grads[i][n] is a VECTOR object containing the gradient on the atom n for the i-th excited state

    #return E_mol, D_mol, E2, sd_basis2, all_grads
    return E_mol, D_mol, S_mol, E2, sd_basis2, all_grads

