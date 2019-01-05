#*********************************************************************************
#* Copyright (C) 2017-2018 Brendan A. Smith, Wei Li, Alexey V. Akimov
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
"""
 This module is designed to convert the results of QE calculations (KS orbital energies and
 time-overlaps in the KS basis) to the generic Hvib matrices, which account for:
 - state reordering;
 - phase corrections;
 - multi-electron wavefunction (Slater determinants) and spin-adaptation
 - scissor operator corrections to energy levels

"""


import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import mapping
import common_utils as comn
import libra_py.tsh as tsh
import libra_py.units as units
import libra_py.hungarian as hungarian



def get_data(params):
    """
    Read the "elementary" overlaps and energies in the basis of KS orbitals

    Required parameter keys:

    params["norbitals"]        [int] - how many lines/columns in the file - the total number of spin-orbitals
    params["active_space"]     [list of ints] - which orbitals we care about (indexing starts with 0)
    params["data_set_path"]    [string] - the path to the directory in which all the files are located
    params["S_re_prefix"]      [string] - prefixes of the files with real part of the MO overlaps at time t
    params["S_re_suffix"]      [string] - suffixes of the files with real part of the MO overlaps at time t
    params["S_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO overlaps at time t
    params["S_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO overlaps at time t
    params["St_re_prefix"]     [string] - prefixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_re_suffix"]     [string] - suffixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_im_prefix"]     [string] - prefixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["St_im_suffix"]     [string] - suffixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["E_re_prefix"]      [string] - prefixes of the files with real part of the MO energies at time t
    params["E_re_suffix"]      [string] - suffixes of the files with real part of the MO energies  at time t
    params["E_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO energies at time t
    params["E_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO energies  at time t
    params["nsteps"]           [int] - how many files to read, starting from index 0
    params["is_pyxaid_format"] [Boolean] - whether the input in the old PYXAID format (just orbital space matrices)

    """

    critical_params = ["norbitals", "active_space", "S_re_prefix", "S_im_prefix",
    "St_re_prefix", "St_im_prefix", "E_re_prefix", "E_im_prefix", "nsteps", "data_set_path" ]

    default_params = { "is_pyxaid_format":False, "S_re_suffix":"_re", "S_im_suffix":"_im",
    "St_re_suffix":"_re", "St_im_suffix":"_im", "E_re_suffix":"_re", "E_im_suffix":"_im"}

    comn.check_input(params, default_params, critical_params)

    norbitals = params["norbitals"]  # the number of orbitals in the input files
    active_space = params["active_space"]
    nstates = len(active_space)
    nsteps = params["nsteps"]

    S, St, E = [], [], [] 

    for i in range(0,nsteps):

        filename_re = params["data_set_path"]+params["S_re_prefix"]+str(i)+params["S_re_suffix"]
        filename_im = params["data_set_path"]+params["S_im_prefix"]+str(i)+params["S_im_suffix"]
        s = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, active_space )
        if params["is_pyxaid_format"]==True:
            s = comn.orbs2spinorbs(s)
        S.append(s)

        filename_re = params["data_set_path"]+params["St_re_prefix"]+str(i)+params["St_re_suffix"]
        filename_im = params["data_set_path"]+params["St_im_prefix"]+str(i)+params["St_im_suffix"]
        st = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, active_space )
        if params["is_pyxaid_format"]==True:
            st = comn.orbs2spinorbs(st)
        St.append(st)

        filename_re = params["data_set_path"]+params["E_re_prefix"]+str(i)+params["E_re_suffix"]
        filename_im = params["data_set_path"]+params["E_im_prefix"]+str(i)+params["E_im_suffix"]
        e = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, active_space )
        if params["is_pyxaid_format"]==True:
            e = comn.orbs2spinorbs(e)
        E.append(e)


    return S, St, E


def get_Lowdin(S):
    """
    Find the S_i_half for the S matrix - alpha and beta components
    """

    nstates = S.num_of_cols/2  # division by 2 because it is a super-matrix

    alp = range(0,nstates)
    bet = range(nstates, 2*nstates)

    S_aa = CMATRIX(nstates, nstates)
    S_bb = CMATRIX(nstates, nstates)

    pop_submatrix(S, S_aa, alp, alp)
    pop_submatrix(S, S_bb, bet, bet)

    is_inv = FullPivLU_rank_invertible(S_aa)
    if is_inv[1] != 1:
        print "Error, S_aa is not invertible, Exiting Program";  sys.exit(0)

    is_inv = FullPivLU_rank_invertible(S_bb)
    if is_inv[1] != 1:
        print "Error, S_bb is not invertible, Exiting Program";  sys.exit(0)

    S_aa_half = CMATRIX(nstates,nstates)
    S_aa_i_half = CMATRIX(nstates,nstates)
    sqrt_matrix(S_aa, S_aa_half, S_aa_i_half)

    S_bb_half = CMATRIX(nstates,nstates)
    S_bb_i_half = CMATRIX(nstates,nstates)
    sqrt_matrix(S_bb, S_bb_half, S_bb_i_half)

    return S_aa_i_half, S_bb_i_half




def apply_normalization(S, St, params):
    """
    Ensure the transition density matrix is computed using the 
    normalized wavefunctions.
    """

    critical_params = [ ]
    default_params = { "do_orthogonalization":0 }
    comn.check_input(params, default_params, critical_params)

    if params["do_orthogonalization"]==1:

        nsteps = len(St)
        nstates = St[0].num_of_cols/2  # division by 2 because it is a super-matrix
        
        alp = range(0,nstates)
        bet = range(nstates, 2*nstates)

        St_aa = CMATRIX(nstates,nstates);  St_ab = CMATRIX(nstates,nstates);
        St_ba = CMATRIX(nstates,nstates);  St_bb = CMATRIX(nstates,nstates);
        
        for i in range(0, nsteps-1):
        
            U1_a, U1_b = get_Lowdin(S[i])    # time n
            U2_a, U2_b = get_Lowdin(S[i+1])  # time n+1          
        
            pop_submatrix(St[i], St_aa, alp, alp);    pop_submatrix(St[i], St_ab, alp, bet)
            pop_submatrix(St[i], St_ba, bet, alp);    pop_submatrix(St[i], St_bb, bet, bet)
        
            St_aa = U1_a.H() * St_aa * U2_a
            St_ab = U1_a.H() * St_ab * U2_b
            St_ba = U1_b.H() * St_ba * U2_a
            St_bb = U1_b.H() * St_bb * U2_b
        
            push_submatrix(St[i], St_aa, alp, alp);   push_submatrix(St[i], St_ab, alp, bet)
            push_submatrix(St[i], St_ba, bet, alp);   push_submatrix(St[i], St_bb, bet, bet)


        
         
def make_cost_mat(orb_mat_inp, en_mat_inp, alpha):
    """
    Makes the cost matrix from a given TDM and information on states' energies

    orb_mat_inp  [CMATRIX(nstates,nstates) or MATRIX(nstates,nstate)] TDM in a given basis
    en_mat_inp   [MATRIX(nstates, nstates)]            Matrix of energies in a given basis [a.u.]
    alpha        [float] Parameter controlling the range of the orbitals that can participate in
                         the reordering. Setting is to 0 makes all orbitals be considered for reordering
                         Setting it to a large number makes the effective number of orbitals participating
                         in the reordering smaller - this can be used to turn off the reordering. [a.u.^-1]

    """

    nstates = orb_mat_inp.num_of_cols
    cost_mat = MATRIX(nstates, nstates)

    for a in xrange(nstates):
        for b in xrange(nstates):

            s = orb_mat_inp.get(a,b)
            s2 = (s*s.conjugate()).real
            dE = (en_mat_inp.get(a,a) - en_mat_inp.get(b,b)).real
            val = s2 * math.exp(-(alpha*dE)**2)
            cost_mat.set(a,b,val)

    return cost_mat



def apply_state_reordering(St, E, params):
    """
    Performs the state's identity reordering in a given basis for all time steps.
    This is reflects in the corresponding changess of the TDM.

    St      [list of CMATRIX] - TDM for each timestep
    E       [list of CMATRIX] - energies of all states at every step
    params  [Python dictionary] - parameters controlling the reordering

    === Required parameter keys: ===

    params["do_state_reordering"]    [int] - option to select an algorithm [default: 2]
                                     Available options:
                                     1 - older version developed by Kosuke Sato, may not the working all the times
                                     2 - Munkres-Kuhn (Hungarian) method

    params["state_reordering_alpha"] [double] - a parameter that controls how many states will be included in the 
                                     reordering

    """

    critical_params = [ ]
    default_params = { "do_state_reordering":2, "state_reordering_alpha":0.0 }
    comn.check_input(params, default_params, critical_params)

    nsteps = len(St)
    nstates = St[0].num_of_cols/2 # division by 2 because it is a super-matrix

    alp = range(0,nstates)
    bet = range(nstates, 2*nstates)

    # Initialize the cumulative permutation as the identity permutation
    perm_cum_aa = intList() # cumulative permutation for alpha spatial orbitals
    perm_cum_bb = intList() # cumulative permutation for beta  spatial orbtials
    for a in xrange(nstates):
        perm_cum_aa.append(a)
        perm_cum_bb.append(a)

    # Current permutation
    perm_t_aa = intList() 
    perm_t_bb = intList()
    for a in xrange(nstates):
        perm_t_aa.append(a)
        perm_t_bb.append(a)


    # Temporary matrices for the Hungarian method
    aa = CMATRIX(nstates, nstates); ab = CMATRIX(nstates, nstates)
    ba = CMATRIX(nstates, nstates); bb = CMATRIX(nstates, nstates)

    en_mat_aa = CMATRIX(nstates, nstates); en_mat_bb = CMATRIX(nstates, nstates)



    for i in range(0, nsteps):

        if params["do_state_reordering"]==1:
            """
            A simple approach based on permuations - but this is not robust
            may have loops
            """
            perm_t = get_reordering(St[i])

            # apply the cumulative permutation  
            update_permutation(perm_t, perm_cum)

            # apply the permutation
            # Because St = <psi(t)|psi(t+dt)> - we permute only columns
            St[i].permute_cols(perm_cum)

            E[i].permute_cols(perm_cum)
            E[i].permute_rows(perm_cum)


        elif params["do_state_reordering"]==2:
            """
            The Hungarian approach
            """

            pop_submatrix(St[i], aa, alp, alp); pop_submatrix(St[i], ab, alp, bet)
            pop_submatrix(St[i], ba, bet, alp); pop_submatrix(St[i], bb, bet, bet)

            # Extract the alpha and beta orbtial energies
            pop_submatrix(E[i], en_mat_aa, alp, alp); pop_submatrix(E[i], en_mat_bb, bet, bet)

            # Permute rows 
            aa.permute_rows(perm_t_aa);   ab.permute_rows(perm_t_aa)
            ba.permute_rows(perm_t_bb);   bb.permute_rows(perm_t_bb)


            # compute the cost matrices for diagonal blocks
            cost_mat_aa = make_cost_mat(aa, en_mat_aa, params["state_reordering_alpha"])
            cost_mat_bb = make_cost_mat(bb, en_mat_bb, params["state_reordering_alpha"])          

            # Solve the optimal assignment problem for diagonal blocks
            res_aa = hungarian.maximize(cost_mat_aa)
            res_bb = hungarian.maximize(cost_mat_bb)
   

            # Convert the list of lists into the permutation object
            for ra in res_aa:
                perm_t_aa[ra[0]] = ra[1]  # for < alpha | alpha > this becomes a new value: perm_t = P_{n+1}
            for rb in res_bb:
                perm_t_bb[rb[0]] = rb[1]  # for < beta | beta > this becomes a new value: perm_t = P_{n+1}   

            # Permute the blocks by col
            aa.permute_cols(perm_t_aa);  ab.permute_cols(perm_t_bb)
            ba.permute_cols(perm_t_aa);  bb.permute_cols(perm_t_bb)


            # Reconstruct St matrix 
            push_submatrix(St[i], aa, alp, alp); push_submatrix(St[i], ab, alp, bet)
            push_submatrix(St[i], ba, bet, alp); push_submatrix(St[i], bb, bet, bet)




def do_phase_corr(cum_phase, St, phase_i):
    """
    This function changes the St matrix according to
    the previous cumulative phases and the current 
    phase correction:

    St -> St = F_n * St * (f_{n+1})^+

    cum_phase [CMATRIX(nstates, 1)]        cumulative phase corrections up to step n (F_n)
    St        [CMATRIX(nstates, nstates)]  input/output TDM to be processed: 
                                           could be alpha-alpha, beta-beta, alpha-beta, 
                                           or beta-alpha sub-blocks
    phase_i   [CMATRIX(nstates, 1)]        the current step phase corrections (f_{n+1})
    """
   
    nstates = St.num_of_rows

    ### Correct the TDM matrix ###
    for a in xrange(nstates):
        for b in xrange(nstates):
            #fab = cum_phase.get(b) * cum_phase.get(b).conjugate() * phase_i.get(b).conjugate()
            fab = cum_phase.get(a) * phase_i.get(b).conjugate()
            St.scale(a,b, fab)



def apply_phase_correction(St, params):
    """
    Perform the phase correction according to:
        
    Akimov, A. V. J. Phys. Chem. Lett, 2018, 9, 6096

    """

    critical_params = [ ]
    default_params = { "do_phase_correction":1 }
    comn.check_input(params, default_params, critical_params)

    nsteps = len(St)
    nstates = St[0].num_of_cols/2  # division by 2 because it is a super-matrix

    alp = range(0,nstates)
    bet = range(nstates, 2*nstates)

    ### Initiate the cumulative phase correction factors ###    
    cum_phase_aa = CMATRIX(nstates,1)  # F(n-1)  cumulative phase
    cum_phase_bb = CMATRIX(nstates,1)  # F(n-1)  cumulative phase
    for a in xrange(nstates):
        cum_phase_aa.set(a, 0, 1.0+0.0j)
        cum_phase_bb.set(a, 0, 1.0+0.0j)


    St_aa = CMATRIX(nstates, nstates); St_ab = CMATRIX(nstates, nstates)
    St_ba = CMATRIX(nstates, nstates); St_bb = CMATRIX(nstates, nstates)

    for i in range(0, nsteps):

        pop_submatrix(St[i], St_aa, alp, alp)
        pop_submatrix(St[i], St_bb, bet, bet)

        if params["do_phase_correction"]==1:

            ### Compute the instantaneous phase correction factors for diag. blocks ###
            phase_i_aa = compute_phase_corrections(St_aa)   # f(i)
            phase_i_bb = compute_phase_corrections(St_bb)   # f(i)       

            ### Do the  phase correstions for the diag. blocks ###
            do_phase_corr(cum_phase_aa, St_aa, phase_i_aa)
            do_phase_corr(cum_phase_bb, St_bb, phase_i_bb)
        
            ### Do the  phase correstions for the off-diag. blocks ###
            do_phase_corr(cum_phase_aa, St_ab, phase_i_bb)
            do_phase_corr(cum_phase_bb, St_ba, phase_i_aa)

            ### Push the corrected diag. blocks to orig. St matrix ###
            push_submatrix(St[i], St_aa, alp, alp);   push_submatrix(St[i], St_ab, alp, bet)
            push_submatrix(St[i], St_ba, bet, alp);   push_submatrix(St[i], St_bb, bet, bet)

            ### Update the cumulative phase correction factors for diag. blocks ###
            for a in xrange(nstates):
                cum_phase_aa.scale(a, 0, phase_i_aa.get(a))
                cum_phase_bb.scale(a, 0, phase_i_bb.get(a))




def sac_matrices(coeff, basis, norbs):

    n_chi = len(coeff)
    n_phi = len(coeff[0])

    P2C = CMATRIX(n_phi, n_chi)
    for j in xrange(n_chi):
        for i in xrange(n_phi):
            P2C.set(i,j,coeff[j][i]*(1.0+0.0j) )

    # Make Sorb
    alp = range(0,norbs/2)
    bet = range(norbs/2, norbs)

    Sorb = CMATRIX(norbs,norbs)
    iden = CMATRIX(norbs/2,norbs/2)
    for i in xrange(norbs/2):
        iden.set(i,i,1.0,0.0)

    push_submatrix(Sorb, iden, alp, alp); push_submatrix(Sorb, iden, alp, bet)
    push_submatrix(Sorb, iden, bet, alp); push_submatrix(Sorb, iden, bet, bet)

    # Compute the overlaps of the SDs:
    #
    Ssd = mapping.ovlp_mat_arb(basis, basis, Sorb)

    # Normalize the Chi wavefunctions #
    norm = (P2C.H() * Ssd * P2C).real()
    for i in xrange(n_chi):
        if norm.get(i,i) > 0.0:
            P2C.scale(-1, i, 1.0/math.sqrt(norm.get(i,i)) )
        else:
            print "Error in CHI normalizaiton: some combination gives zero norm\n"
            sys.exit(0)


    return P2C


def scale_H_vib(hvib, en_gap, dNAC, sc_nac_method):
    """
    This function scales the energies and NACs in the vibrionic Hamiltonian
    in the Chi basis.   
    hvib   [list of CMATRIX objects] = CMATRIXlist of vibronic hamiltonians in the Chi basis
    en_gap [float]  = The desired energy gap (E_1 - E_0), for the Chi basis
    dNAC   [list of lists of (list, float)] = The scaling terms by which specific nacs will 
                                              be scaled datatype = list of lists of (list, float)

                                                         [  [ [i,j], val ], ...  ]          

                                             n and n+1 are the col (and thereby row) indicies of 
                                             the nacs to be scaled by the value val 
    sc_nac_method [int] = The method used to scale NACs in the Chi basis, chosen by the user.
                          If sc_nac_method = 1, then the NACs are scaled by the ivnerse of the
                          magnitude of the change in energy, according to Lin et al.

                          Reference:
                          Lin, Y. & Akimov, A. V. J. Phys. Chem. A (2016) 
    """

    traj_len = len(hvib)

    # Do the scaling
    nrows = hvib[0].num_of_rows
    ncols = hvib[0].num_of_cols   
    for i in xrange(traj_len):

        prev_gap = hvib[i].get(1,1) - hvib[i].get(0,0)
        shift = en_gap - prev_gap 

        # Scale the energy gap, by adding the scaling factors to all excited states
        # This is to keep the ground state enegy equal to 0.0 by definition
        for j in xrange(1,ncols):
            hvib[i].add(j, j, shift)

        if sc_nac_method == 0:

            # Scales nacs manually as set by user
            for it in dNAC:
                a,b = it[0][0], it[0][1]
                val = it[1]
                hvib[i].scale(a,b, val)
                hvib[i].scale(b,a, val)

        elif sc_nac_method == 1:
            # Scales nacs by the inverse of the change in energy gap
            for j in xrange(1,ncols):
                scl_nac = prev_gap / (shift + prev_gap)
                hvib[i].scale(0, j, scl_nac)
                hvib[i].scale(j, 0, scl_nac)

    return hvib                    


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




def run(params):
    """
    The procedure to converts the results of QE calculations (KS orbital energies and
    time-overlaps in the KS basis) to the generic Hvib matrices, which account for:   
    - state reordering;
    - phase corrections;
    - multi-electron wavefunction (Slater determinants) and spin-adaptation
    - scissor operator corrections to energy levels

    === Required parameter keys: ===

    params["Phi_basis"]        [list of lists of ints] - define the Slater Determinants basis
    params["P2C"]              [list of lists of complex] - define the superpositions to SDs to get spin-adapted functions
    params["Phi_dE"]           [list of doubles] - define corrections of the SAC state energies
    params["dt"]               [double] - nuclear dynamics integration time step [in a.u. of time, default: 41.0]

    params["do_state_reordering"]     [int] - option to do the state reordering [default: 2]: 
                                 0 - no state reordering
                                 1 - older method (may be wrong)
                                 2 - Hungarian algorithm

    params["state_reordering_alpha"]  [double] - the energy window parameter in the Hungarian approach  [default: 0.00]
    params["do_phase_correction"]     [int] - option to do the phase correction [default: 1]
                                 0 - don't do 
                                 1 - do


    === Required by the get_data() ===

    params["norbitals"]        [int] - how many lines/columns in the file - the total number of spin-orbitals
    params["active_space"]     [list of ints] - which orbitals we care about (indexing starts with 0)
    params["nsteps"]           [int] - how many files to read, starting from index 0
    params["is_pyxaid_format"] [Boolean] - whether the input in the old PYXAID format (just orbital space matrices)   
    params["data_set_paths"]   [string] - where the input files are located
    params["S_re_prefix"]      [string] - prefixes of the files with real part of the MO overlaps at time t
    params["S_re_suffix"]      [string] - suffixes of the files with real part of the MO overlaps at time t
    params["S_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO overlaps at time t
    params["S_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO overlaps at time t
    params["St_re_prefix"]     [string] - prefixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_re_suffix"]     [string] - suffixes of the files with real part of the MO overlaps at times t and t' = t+dt
    params["St_im_prefix"]     [string] - prefixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["St_im_suffix"]     [string] - suffixes of the files with imaginary part of the MO overlaps at times t and t' = t+dt
    params["E_re_prefix"]      [string] - prefixes of the files with real part of the MO energies at time t
    params["E_re_suffix"]      [string] - suffixes of the files with real part of the MO energies  at time t
    params["E_im_prefix"]      [string] - prefixes of the files with imaginary part of the MO energies at time t
    params["E_im_suffix"]      [string] - suffixes of the files with imaginary part of the MO energies  at time t


    === Define the output ===
    params["output_set_paths"] [string] - where the resulting Hvib files are to be stored [default: the same as the input paths]
    params["Hvib_re_prefix"]   [string] - prefixes of the output files with real part of the vibronic Hamiltonian at time t
    params["Hvib_re_suffix"]   [string] - suffixes of the output files with real part of the vibronic Hamiltonian at time t
    params["Hvib_im_prefix"]   [string] - prefixes of the output files with imaginary part of the vibronic Hamiltonian at time t
    params["Hvib_im_suffix"]   [string] - suffixes of the output files with imaginary part of the vibronic Hamiltonian at time t
     
    """

    critical_params = [ "P2C", "Phi_basis", "Phi_dE", "data_set_paths" ]
    default_params = { "dt":41.3413,  "output_set_paths":params["data_set_paths"], 
                       "Hvib_re_prefix":"Hvib_", "Hvib_im_prefix":"Hvib_",
                       "Hvib_re_suffix":"_re", "Hvib_im_suffix":"_im",
                       "do_state_reordering":2, "state_reordering_alpha":0.0,
                       "do_phase_correction":1
                     }
    comn.check_input(params, default_params, critical_params)
 
    if(len(params["data_set_paths"]) != len(params["output_set_paths"])):
        print "Error: Input and output sets paths should have equal number of entries\n"
        print "len(params[\"data_set_paths\"]) = ", len(params["data_set_paths"])
        print "len(params[\"output_set_paths\"]) = ", len(params["output_set_paths"])
        print "Exiting...\n"
        sys.exit(0)

    dt = params["dt"]

    """
    1. Read in the "elementary" overlaps and energies in the basis of KS orbitals
    2. Apply state reordering to KS
    3. Apply phase correction to KS
    4. Construct the Hvib in the basis of Slater determinants
    5. Convert the Hvib to the basis of symmery-adapted configurations (SAC)
    """

    P2C = sac_matrices(params["P2C"], params["Phi_basis"], params["norbitals"])
    H_vib = []
    ndata = len(params["data_set_paths"])

    for idata in xrange(ndata):

        prms = dict(params)
        prms.update({"data_set_path":params["data_set_paths"][idata]})    
        S_dia_ks, St_dia_ks, E_dia_ks = get_data(prms)  

        apply_normalization(S_dia_ks, St_dia_ks, prms)
        apply_state_reordering(St_dia_ks, E_dia_ks, prms)    
        apply_phase_correction(St_dia_ks, prms)
        
        nsteps = len(St_dia_ks)
        nstates = len(prms["P2C"])
        
        Hvib = []
        for i in xrange(nsteps):
        
            # Hvib in the basis of SDs
            hvib = compute_Hvib(prms["Phi_basis"], St_dia_ks[i], E_dia_ks[i], prms["Phi_dE"], dt) 
        
            # SAC
            Hvib.append( P2C.H() * hvib * P2C )

        #========== Scale H_vibs ===============
        # Document the params keys - then we may uncomment this
        #if params["do_scale"] == 1:
        #    scale_H_vib(Hvib, params["Chi_en_gap"], params["NAC_dE"], params["sc_nac_method"])

        
        # Output the resulting Hamiltonians
        for i in xrange(nsteps):

            re_filename = prms["output_set_paths"][idata] + prms["Hvib_re_prefix"] + str(i) + prms["Hvib_re_suffix"]
            im_filename = prms["output_set_paths"][idata] + prms["Hvib_im_prefix"] + str(i) + prms["Hvib_im_suffix"]        

            Hvib[i].real().show_matrix(re_filename)
            Hvib[i].imag().show_matrix(im_filename)

            # Make time-derivative overlap matriices for chi basis and print
            #St_phi = mapping.ovlp_mat_arb(params["Phi_basis"], params["Phi_basis"], St_dia_ks[i])
            #St_chi = P2C.H() * St_phi * P2C

            #re_filename = prms["output_set_paths"][idata] + prms["St_SD_re_prefix"] + str(i) + prms["St_SD_re_suffix"]
            #St_chi.real().show_matrix(re_filename)

        H_vib.append(Hvib)        
        

    return H_vib

