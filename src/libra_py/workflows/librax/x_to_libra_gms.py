#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#* 
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file x_to_libra_gms.py 
# This module implements the functions that extract parameters from the gamess output file:
# atomic forces , molecular energies, molecular orbitals, and atomic basis information.
# The forces are used for simulating Classical MD on Libra 
# and the others for calculating time-averaged energies and Non-Adiabatic Couplings(NACs).

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *
from extract_gms import *
#from overlap import *
#from hamiltonian_el import *
from moment import *
from misc import *
from spin_indx import *
import reorder_matrices
import unavoided_tmp # set temporally

def exe_gamess(params):
    ##
    # This is a function that call GAMESS execution on the compute node
    # \param[in] params Input data containing all manual settings and some extracted data.
    #
    # Used in main.py/main and md.py/run_MD

    inp = params["gms_inp"]
    out = params["gms_out"]
    nproc = params["nproc"]

    scr_dir = params["scr_dir"]
    rungms = params["rungms"]
    VERNO = params["VERNO"]

    # set environmental variables for GAMESS execution
    os.environ["SCR"] = scr_dir
    os.environ["USERSCR"] = scr_dir
    os.environ["GMSPATH"] = params["GMSPATH"]

    # create scratch directory
    os.system("mkdir %s" % (scr_dir))

    #os.system("/usr/bin/time rungms.slurm %s 01 %s > %s" % (inp,nproc,out))
    os.system("/usr/bin/time %s %s %s %s > %s" % (rungms,inp,VERNO,nproc,out))

    # delete the files except input and output ones to do another GAMESS calculation.
    os.system("rm *.dat")
    os.system("rm -r %s" % (scr_dir))


def gamess_to_libra(params, ao, E, sd_basis, active_space,suff):
    ## 
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params         contains input parameters , in the dictionary form
    # \param[in,out] ao         atomic orbital basis at "t" old
    # \param[in,out] E          total excitation energies at "t" old
    # \param[in] sd_basis       Basis of Slater determinants at "t" old (list of CMATRIX object).
    #                           In the present implementation, it contains a single determinant
    # \param[in] active_space   A list of indices (starting from 1) of the MOs to include in calculations (and to read from the QE output files)  
    # \param[in] suff           A suffix to add to the name of the output files; this suffix is now considered to be of a string type 
    #                           - so you can actually encode both the iteration number (MD timestep),
    #                           the nuclear cofiguration (e.g. trajectory), and any other related information
    # Returned datas are explained above the return line.
    #
    # Used in: md.py/run_MD

    flag_ao = params["flag_ao"]
    sd_basis2 = SDList()    # this is a list of SD objects. Eeach represents a Slater Determinant
    nstates = len(params["excitations"])
    sz = len(active_space)

    rnd = Random()

    # 2-nd file - time "t+dt"  new
    label, Q, R, Grad, E2, c2, ao2, nel = gms_extract(params["gms_out"],params["excitations"],params["min_shift"],active_space,params["debug_gms_unpack"])

    #e2 = MATRIX(E2)
    homo = params["nel"]/2 +  params["nel"] % 2

    for ex_st in xrange(nstates): 
        mo_pool_alp = CMATRIX(c2)
        mo_pool_bet = CMATRIX(c2)
        alp,bet = index_spin(params["excitations"][ex_st],active_space, homo)

        # use excitation object to create proper SD object for different excited state
        sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet))
        sd_basis2.append(sd)

    # Gradients
    # in this implementation (CPA), the gradients on all excited states are the same

    all_grads = []
    for i in xrange(nstates):
        grd = []
        for g in Grad:
            grd.append(VECTOR(g))
        all_grads.append(grd)
    
    #t = Timer()

    ##P11, P22, P12, P21 = overlap(ao,ao2,sd_basis[0],sd_basis2,params["basis_option"])
    # calculate overlap matrix of Slater determinant basis states
    P11 = SD_overlap(sd_basis,  sd_basis)
    P22 = SD_overlap(sd_basis2, sd_basis2)
    P12 = SD_overlap(sd_basis,  sd_basis2)
    P21 = SD_overlap(sd_basis2, sd_basis)
    #print "Time to compute in SD_overlap= ",t.show(),"sec"

    p0 = range(nstates)
    #perm = unavoided.get_reordering(P12)
    perm = unavoided_tmp.get_reordering(P12) # temporally
    if p0 != perm:
        print "trivial crossings occured!"
        print "perm is", perm
        print "P12 is"; P12.show_matrix()

        reorder_matrices.reorder(perm,P12,E2)
        P21 = P12.H()

    # calculate transition dipole moment matrices in the MO basis:
    # mu_x = <i|x|j>, mu_y = <i|y|j>, mu_z = <i|z|j>
    # this is done for the "current" state only    
    # ********* we use here only SD basis, which means we will calculate dipole moments
    # ********* based on SD basis.

    mu_x, mu_y, mu_z = CMATRIX(sz,sz),CMATRIX(sz,sz),CMATRIX(sz,sz)
    mu = [mu_x, mu_y, mu_z] # initialize mu

    # *************************************************************************
    # Modify here if calculation on dipole moments with SD_overlap can be done.
    # *************************************************************************
    #if flag_ao == 1:
        #t.start()
    #    mu_x, mu_y, mu_z = transition_dipole_moments(ao2,sd_basis2)
    #    mu = [mu_x, mu_y, mu_z] # now mu is defined as a CMATRIX list.
        #t.stop()
        #print "Time to compute in dipole moment= ",t.show(),"sec"
    # *************************************************************************

    if params["debug_mu_output"]==1:
        print "mu_x:";    mu_x.show_matrix()
        print "mu_y:";    mu_y.show_matrix()
        print "mu_z:";    mu_z.show_matrix()
 
    if params["debug_densmat_output"]==1:
        print "P11 and P22 matrices should show orthogonality"
        print "P11 is";    P11.show_matrix()
        print "P22 is";    P22.show_matrix()
        print "P12 and P21 matrices show overlaps of SDs for different molecular geometries "
        print "P12 is";    P12.show_matrix()
        print "P21 is";    P21.show_matrix()

    # Here, explicit computation would be more convinient than using outer functions (in hamiltonian_el.py).
    E_ave = 0.50 * ( E + E2 )
    nac = 0.50/params["dt_nucl"] * ( P12 - P21 )

    # shift exciting energies according to "shift_E" list.
    if "shift_E" in params:
        eV_to_au = 0.036749
        #print "before shifting E(0,0) is %f " % E_ave.get(0,0)
        for i in xrange(nstates):
            etmp = E_ave.get(i,i) + params["shift_E"][i]*eV_to_au
            E_ave.set(i, i, etmp)
        #print "after shifting E(0,0) is %f " % E_ave.get(0,0)

    # here, E_ave and nac are printed for debugging 
    #if 0==1:
    #    E_ave.show_matrix(params["mo_ham"] + "E_ave_" + suff)
    #    nac.real().show_matrix(params["mo_ham"] + "nac_real_" + suff)
    #    nac.imag().show_matrix(params["mo_ham"] + "nac_imag_" + suff)
    #E_mol_red.show_matrix(params["mo_ham"] + "reduced_re_Ham_" + suff)
    #D_mol.show_matrix(params["mo_ham"] + "reduced_im_Ham_" + suff)
    # ********** "CMATRIX.show_matrix(filename)" is not exported ****** 

    # store "t+dt"(new) parameters on "t"(old) ones
    for i in range(0,len(ao2)):
        ao[i] = AO(ao2[i])
    E = MATRIX(E2)  # at time t+dt

    # useless lines: nac is already defined as CMATRIX.
    #nac = CMATRIX(D_mol.num_of_rows, D_mol.num_of_cols)
    #for i in xrange(D_mol.num_of_rows):
    #    for j in xrange(D_mol.num_of_cols):
    #        nac.set(i,j,D_mol.get(i,j),0.0)

    # Returned data:
    # E_ave : the matrix of the total excitation energy averaged over energies at "t" and "t+dt"
    # nac (CMATRIX): the matrix of the NACs computed with SD orbitals. Same dimension as E_ave
    # sd_basis2 : the list of reduced SD (active space orbitals), representing all computed states
    # all_grads: all_grads[i][k] - the gradient w.r.t. to k-th nucleus of i-th excitation state
    # mu : mu[i] transition dipole moment of i-th DOF. (mu_x, mu_y, mu_z)

    return E, E_ave, nac, sd_basis2, all_grads, mu

