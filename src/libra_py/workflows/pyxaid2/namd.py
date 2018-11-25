#*********************************************************************************
#* Copyright (C) 2017 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  
#
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

import mapping
import libra_py.workflows.common_utils as comn


def show_matrix_splot(X, filename):
    ncol, nrow = X.num_of_cols, X.num_of_rows

    line = ""
    for i in xrange(nrow):
        for j in xrange(ncol):
            val = X.get(i,j)
            if i==j:
                val = 0.0
            line = line + "%4i %4i %8.5f \n" % (i, j, val)
        line = line + "\n"

    f = open(filename, "w")
    f.write(line)
    f.close()
 
    

def add_printout(i, pop, filename):
    # pop - CMATRIX(nstates, 1)

    f = open(filename,"a")
    line = "step= %4i " % i    

    tot_pop = 0.0
    for st in xrange(pop.num_of_cols):
        pop_o = pop.get(st,st).real
        tot_pop = tot_pop + pop_o
        line = line + " P(%4i)= %8.5f " % (st, pop_o)
    line = line + " Total= %8.5f \n" % (tot_pop)
    f.write(line)
    f.close()


def compute_Hvib(basis, St_ks, E_ks, dE, dt):
    """
    Basis - list of list of lists of integers
    St_ks - time overlap (CMATRIX) of the KS spin-orbitals
    E_ks - energies of KS spin-orbitals at the mid-point
    dE - energy corrections to the SD orbitals ("scissor operator")
    dt - the timestep for MD integrations

    Returns: The Vibronic Hamiltonian
    """

    St    = mapping.ovlp_mat_arb(basis, basis, St_ks) 
    H_el  = mapping.energy_mat_arb(basis, E_ks, dE)
    H_vib = H_el - (0.5j/dt)*(St-St.H())

    return H_vib


def run_namd(params):

    rnd = Random()
    T = params["T"]
    kT = units.kB * T
    init_time = params["init_time"]  # integer
    dt = params["dt"]
    len_traj = params["len_traj"]
    num_sh_traj = params["num_sh_traj"]
    sh_method = params["sh_method"]
    do_collapse = params["do_collapse"]
    bolt_opt = params["Boltz_opt"]

    psi_dia_ks = params["psi_dia_ks"];   nst_dia_ks = 2*len(psi_dia_ks)
    active_space = range(0, nst_dia_ks) #params["active_space"]

    Phi_basis = params["Phi_basis"]
    P2C = params["P2C"];                 nst_dia_sac = P2C.num_of_cols
  
    # ------------------read and store the projection and energies------------------
    H_vib = []
    cum_phase = None  # F(n-1)

    for i in range(0,len_traj):

        ##############################################################################
        # Read in the "elementary" overlaps and energies - in the basis of KS orbitals
        ##############################################################################       

        filename_re = params["E_dia_ks_re_prefix"]+str(i)+params["E_dia_ks_re_suffix"]
        filename_im = params["E_dia_ks_im_prefix"]+str(i)+params["E_dia_ks_im_suffix"]
        E_dia_ks = comn.get_matrix(nst_dia_ks, nst_dia_ks, filename_re, filename_im, active_space )

        filename_re = params["S_dia_ks_re_prefix"]+str(i)+params["S_dia_ks_re_suffix"]
        filename_im = params["S_dia_ks_im_prefix"]+str(i)+params["S_dia_ks_im_suffix"]
        S_dia_ks = comn.get_matrix(nst_dia_ks, nst_dia_ks, filename_re, filename_im, active_space )

        filename_re = params["St_dia_ks_re_prefix"]+str(i)+params["St_dia_ks_re_suffix"]
        filename_im = params["St_dia_ks_im_prefix"]+str(i)+params["St_dia_ks_im_suffix"]
        St_dia_ks = comn.get_matrix(nst_dia_ks, nst_dia_ks, filename_re, filename_im, active_space )

        sz = St_dia_ks.num_of_rows 


        ### Perform state reordering (must be done before the phase correction) ###

        if params["do_state_reordering"]==1:
            """
            A simple approach based on permuations - but this is not robust
            may have loops
            """

            perm = get_reordering(St_dia_ks)

            St_dia_ks.permute_cols(perm)
            St_dia_ks.permute_rows(perm)

            E_dia_ks.permute_cols(perm)
            E_dia_ks.permute_rows(perm)

        elif params["do_state_reordering"]==2:
            """
            The Hungarian approach
            """
            # construct the cost matrix
            alp = 0.0
            if "state_reordering_alpha" in params.keys():
                alp = params["state_reordering_alpha"]

            cost_mat = MATRIX(sz, sz)
            for a in xrange(sz):
                for b in xrange(sz):
                    s = St_dia_ks.get(a,b)
                    s2 =  (s*s.conjugate()).real
                    dE = (E_dia_ks.get(a,a)-E_dia_ks.get(b,b)).real/kT
                    val = s2 * math.exp(-alp*dE**2)
                    cost_mat.set(a,b, val)

            # run the assignment calculations
            res = hungarian.maximize(cost_mat)

            # convert the list of lists into the permutation object
            perm = intList() 
            for a in xrange(sz):
                perm.append(a)
            for r in res:
                perm[r[0]] = r[1]

            # apply the permutation
            St_dia_ks.permute_cols(perm)
            St_dia_ks.permute_rows(perm)

            E_dia_ks.permute_cols(perm)
            E_dia_ks.permute_rows(perm)


        ### Perform phase correction ###
        if params["do_phase_correction"]==1:
            ### Initiate the cumulative phase correction factors ###
            if i==0:
                cum_phase = CMATRIX(sz,1)
                for a in xrange(sz):
                    cum_phase.set(a, 0, 1.0+0.0j)

            ### Compute the instantaneous phase correction factors ###
            phase_i = compute_phase_corrections(St_dia_ks)   # f(i)

            ### Correct the overlap matrix ###
            for a in xrange(sz):   
                for b in xrange(sz): 
                    fab = cum_phase.get(b) * cum_phase.get(b).conjugate() * phase_i.get(b).conjugate()
                    St_dia_ks.scale(a,b, fab)

            ### Update the cumulative phase correction factors ###
            for a in xrange(sz):
                cum_phase.scale(a, 0, phase_i.get(a))
          

        ### Done with the phase correction ###


        # Printing what we just extracted for t = 0
        #"""
        if i == 0:
            print "\nE_dia_ks = "
            E_dia_ks.real().show_matrix()
            print "\nS_dia_ks = "
            S_dia_ks.real().show_matrix()  
            print "\nSt_dia_ks = "
            St_dia_ks.real().show_matrix()
            Hvib = compute_Hvib(Phi_basis, St_dia_ks, E_dia_ks, params["Phi_dE"], params["dt"] )
            print "\nHvib_Phi_dia.imag()\n"
            Hvib.imag().show_matrix()
            print "\nHvib_Phi_dia.real()\n"
            Hvib.real().show_matrix()
            print "\n Making H_vib (Chi_basis)"
            Hvib = P2C.H() * Hvib * P2C
            print "\nHvib_Chi_dia.imag()\n"
            Hvib.imag().show_matrix()
            print "\nHvib_Chi_dia.real()\n"
            Hvib.real().show_matrix()
            #sys.exit(0)        
        #"""

        H_vib.append( compute_Hvib(Phi_basis, St_dia_ks, E_dia_ks, params["Phi_dE"], params["dt"]) )
        H_vib[i] = P2C.H() * H_vib[i] * P2C

    #========== Initialize the wavefunction amplitudes ===============
    init_Chi = params["init_Chi"]    

    # TD-SE coefficients
    Coeff_Chi = [];  
    Coeff_Phi = [];

    for tr in xrange(num_sh_traj):

        Coeff_Chi.append(CMATRIX(nst_dia_sac, 1)); 
        Coeff_Chi[tr].set(init_Chi, 1.0, 0.0)
        Coeff_Phi.append(P2C*Coeff_Chi[tr]) 

        if tr==0:
            print "Coeff_Chi = "; Coeff_Chi[tr].show_matrix()
            print "Coeff_Phi = "; Coeff_Phi[tr].show_matrix()


    #========== Initialize the SE density matrices ===============
    # Density matrices
    denmat_Chi_se = tsh.amplitudes2denmat(Coeff_Chi)
    denmat_Phi_se = tsh.amplitudes2denmat(Coeff_Phi)

    print "denmat_Chi_se = "; denmat_Chi_se[0].show_matrix()
    print "denmat_Phi_se = "; denmat_Phi_se[0].show_matrix()

    #========== Initialize the state for adiabatic sh ===============
    istate  = []  # in the Chi basis (referring to spin-adaapted states)
    for tr in xrange(num_sh_traj):
        prob = tsh.denmat2prob(denmat_Chi_se[tr])
        st = tsh.set_random_state(prob, rnd.uniform(0.0, 1.0)) 
        istate.append(st)

    #========== Initialize output files ===============
    f = open("_pop_Chi_se.txt","w");   f.close()
    f = open("_pop_Chi_sh.txt","w");   f.close()
    f = open("_pop_Phi_se.txt","w");   f.close()
    f = open("_pop_Phi_sh.txt","w");   f.close()


    #=============== Entering Dynamics !! =========================  
    #=============== Now handle the electronic structure ==========
    for i in xrange(init_time, init_time + len_traj):

        #============== TD-SE and surface hopping =================
        for tr in xrange(num_sh_traj):

            # Coherent evolution for Chi
            propagate_electronic(dt, Coeff_Chi[tr], H_vib[i])  # propagate in the diabatic basis (Chi)

            ksi  = rnd.uniform(0.0, 1.0);
            ksi2 = rnd.uniform(0.0, 1.0)

            # Surface hopping in Chi basis
            istate[tr] = tsh.hopping(Coeff_Chi[tr], H_vib[i], istate[tr], sh_method, do_collapse, ksi, ksi2, dt, T, bolt_opt)
            Coeff_Phi[tr] = P2C*Coeff_Chi[tr]

        #============== Analysis of Dynamics  =====================
        # Update SE-derived density matrices
        denmat_Chi_se = tsh.amplitudes2denmat(Coeff_Chi)
        denmat_Phi_se = tsh.amplitudes2denmat(Coeff_Phi)

        # Use SE-derived density matrices to make SH-derived density matrices
        denmat_Chi_sh = []
        denmat_Phi_sh = []
        for tr in xrange(num_sh_traj):

            denmat_Chi_sh.append(CMATRIX(nst_dia_sac, nst_dia_sac))
            denmat_Chi_sh[tr].set(istate[tr],istate[tr], 1.0, 0.0)
            denmat_Phi_sh.append( P2C*denmat_Chi_sh[tr]*P2C.H() )

        ave_pop_Chi_sh, ave_pop_Chi_se = tsh.ave_pop(denmat_Chi_sh, denmat_Chi_se)
        ave_pop_Phi_sh, ave_pop_Phi_se = tsh.ave_pop(denmat_Phi_sh, denmat_Phi_se)

        #===================== Print out ==============================
        add_printout(i, ave_pop_Chi_sh, "_pop_Chi_sh.txt")
        add_printout(i, ave_pop_Chi_se, "_pop_Chi_se.txt")

        add_printout(i, ave_pop_Phi_sh, "_pop_Phi_sh.txt")
        add_printout(i, ave_pop_Phi_se, "_pop_Phi_se.txt")


