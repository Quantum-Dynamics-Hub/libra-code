#*********************************************************************************
#* Copyright (C) 2017-2019 Brendan A. Smith, Wei Li, Alexey V. Akimov
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
.. module:: step3
   :platform: Unix, Windows
   :synopsis: This module is designed to convert the results of QE calculations 
       (KS orbital energies and time-overlaps in the KS basis) to the generic Hvib
       matrices, which account for:

           - state reordering;
           - phase corrections;
           - multi-electron wavefunction (Slater determinants) and spin-adaptation
           - scissor operator corrections to energy levels

.. moduleauthor:: Brendan A. Smith, Wei Li, and Alexey V. Akimov

"""


import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import mapping
import util.libutil as comn
import libra_py.tsh as tsh
import libra_py.units as units
import libra_py.data_read as data_read
import libra_py.hungarian as hungarian


def get_step2_data(_params):
    """
    A light function to obtain the step2 data

    Args:
        params ( dictionary ): Control paramerter of this type of simulation. Can include the follwing keys:

            * **params["basis"]** ( string ): describes if one is using either the spin-diabatic (non spin-orbit coupling) 
                                              or is using the spin-adiabatic (spin orbit-coupling)
    """

    params = dict(_params)
    spin_basis = params["basis"]

    # Fetching the overlap matricies
    params.update({ "data_re_prefix" : "S_"+spin_basis+"_ks_", "data_re_suffix" : "_re",
                    "data_im_prefix" : "S_"+spin_basis+"_ks_", "data_im_suffix" : "_im"  } )
    S = data_read.get_data_sets(params)
    #print ("\n")
    #S[0][0].show_matrix()
    #print (S[0][0].get(0,0))
    #sys.exit(0)

    # Fetching the time-derivative overlap matricies
    params.update({ "data_re_prefix" : "St_"+spin_basis+"_ks_", "data_re_suffix" : "_re",
                    "data_im_prefix" : "St_"+spin_basis+"_ks_", "data_im_suffix" : "_im"  } ) 
    St = data_read.get_data_sets(params)
    #print ("\n")
    #St[0][0].show_matrix()
    #sys.exit(0)
 
    # Fetching the vibronic Hamiltonian matricies
    params.update({ "data_re_prefix" : "hvib_"+spin_basis+"_", "data_re_suffix" : "_re",
                    "data_im_prefix" : "hvib_"+spin_basis+"_", "data_im_suffix" : "_im"  } ) 
    Hvib_ks = data_read.get_data_sets(params)
    #print ("\n")
    #Hvib_ks[0][0].show_matrix()
    #sys.exit(0)

    return S, St, Hvib_ks




def print_SD_basis(SD_basis):
    """
    Just a light function to print the SD basis

    Args:
        SD_basis ( list of lists of ints ): Slater determinant basis in terms of Kohn-Sham orbital indicies

            Possible ground state configurations

            Ex. 1. SD_basis[0] = [ 5, -15 ]
            Ex. 2. SD_basis[0] = [ 2, 3, 4, 5, -12, -13, -14, -15 ]

            Possible ground state configurations

            Ex. 1. SD_basis[N > 0] = [ 6, -15 ]
            Ex. 2. SD_basis[N > 0] = [ 2, 3, 4, 6, -12, -13, -14, -15 ]

            Where we have 10 alpha and beta spin orbitals, and the cbm index for the alpha spin-channel is 5
            and the cbm index for the beta spin-channel is 15. 
    """

    for i in range(len(SD_basis)):
        if i == 0:
            print (" GS: ", SD_basis[i])
        else:
            print (" ES "+str(i)+": ", SD_basis[i]) 
    print ("WARNING: The Slater determinant bases above may not be sorted based on their energies\n")




def sort_SD_energies(Hvib):
    """
    This function goes into the Hvib (SD basis) files and sorts the energies
    For each Hvib file, we are going to obtain a list of lists of orbital index and energy pair
    These orbital index and energy pairs for each step will be sorted based on energies  

    Args:
        Hvib ( list of lists of CMATRIX objects ): vibronic Hamiltonian in the Slater determinant basis

    """

    ntraj  = len(Hvib)
    nsnaps = len(Hvib[0])   
    nSD = Hvib[0][0].num_of_cols

    orbital_index_energy_pairs = []

    for traj in range(ntraj):

        orbital_index_energy_pairs.append( [] )

        for snap in range(nsnaps):

            index_energy_pairs = []

            for sd_index in range(nSD):
                index_energy_pairs.append( [ sd_index, Hvib[ traj ][ snap ].get( sd_index, sd_index ).real ]  )

            sorted_index_energy_pairs = merge_sort( index_energy_pairs  ) 
            orbital_index_energy_pairs[traj].append( sorted_index_energy_pairs )

    return orbital_index_energy_pairs




def output_sorted_Hvibs(Hvib, orbital_index_energy_pairs):
    """
    This function outputs the vibronic Hamiltonians in the SD basis according to their sorted order

    Args:
        Hvib ( list of lists of CMATRIX objects ): vibronic Hamiltonian in the Slater determinant basis
        orbital_index_energy_pairs ( list of lists of lists of lists ): orbital index and energy pair lists for each step and nuclear trajectory
                               Ex. orbital_index_energy_pairs[i][j][k][0] = kth slater determinant index at the jth step on the ith nuclear trajectory   
                               Ex. orbital_index_energy_pairs[i][j][k][1] = kth slater determinant energy at the jth step on the ith nuclear trajectory   

    """

    ntraj  = len(orbital_index_energy_pairs)
    nsnaps = len(orbital_index_energy_pairs[0])
    nSD    = len(orbital_index_energy_pairs[0][0])

    Hvibs_sorted   = []

    for traj in range( ntraj ):

        rd_sorted = "res_sorted_traj"+str(traj)+""
        os.system("rm -r "+rd_sorted)
        os.system("mkdir "+rd_sorted)

        Hvibs_sorted.append( [] )

        for snap in range( nsnaps ):
            Hvib_sorted  = CMATRIX( nSD, nSD)
            for i in range( nSD ):
                for j in range( nSD ):
                    a = orbital_index_energy_pairs[ traj ][ snap ][ i ][ 0 ]
                    b = orbital_index_energy_pairs[ traj ][ snap ][ j ][ 0 ]
                    Hvib_sorted.set( i, j, Hvib[ traj ][ snap ].get( a, b ) )

            Hvibs_sorted[ traj ].append( Hvib_sorted )
            Hvibs_sorted[ traj ][ snap ].real().show_matrix("%s/Hvib_sorted_%i_re" % (rd_sorted, snap))
            Hvibs_sorted[ traj ][ snap ].imag().show_matrix("%s/Hvib_sorted_%i_im" % (rd_sorted, snap))

    return Hvibs_sorted




def build_SD_basis(data_dim, cbm_alpha_index, alpha_include, cbm_beta_index, beta_include, excitation_type):
    """
    Builds a Slater Determinant basis based on the indexing notation scheme used in Libra

    Args:
        data_dim (int): how many rows or columns in the vibronic Hamiltonian matrix. This will be an even number becuase the 
                        number of alpha orbtials should equal the number of beta orbitals
        cbm_(alpha/beta)_index (int): index of VBM (or HOMO) in the matrix of the vibronic Hamiltonian 
                                      (row or column index). Note, this index is from 1 
        (alpha/beta)_include (int): how many orbitals to include from the cbm_(alpha/beta)_index
        excitation_type (int):  0: Make SDs with beta electrons excited
                                1: Make SDs with alpha electrons excited
                                2: Make two sets of SDs, one for beta and one for alpha electrons excited  
    """

    num_same_spin_orbitals = int( data_dim / 2 )

    alpha_electrons = []
    beta_electrons  = []

    # Check for potential errors
    if data_dim % 2 > 0:
        print ("The dimensions of your vibronic Hamiltonian matrix (or whatever data you wish to pass to the function build_SD) must be even")
        print ("Exiting now ..")
        sys.exit(0)

    if cbm_alpha_index > num_same_spin_orbitals:
        print ("You must have the same number of alpha and beta spin-orbitals in your basis. The index of the CBM for the alpha spin-orbitals\
 cannot be greater than 1/2 the total number of spin orbitals ")
        print ("Exiting now")
        sys.exit(0)

    elif cbm_alpha_index <= 0:
        print ("The index of the CBM for the alpha spin-orbitals must be > 0")
        print ("Exiting now")
        sys.exit(0)

    elif cbm_beta_index <= num_same_spin_orbitals:
        print ("You must have the same number of alpha and beta spin-orbitals in your basis. The index of the CBM for the beta spin-orbitals\
 cannot be less than 1/2 the total number of spin orbitals ")
        print ("Exiting now")
        sys.exit(0)

    elif cbm_beta_index > 2*num_same_spin_orbitals:
        print ("The index of the CBM for the beta spin-orbitals must be < total number of spin-orbitals")
        print ("Exiting now")
        sys.exit(0)

    if cbm_alpha_index + alpha_include + 2 > num_same_spin_orbitals:
        print ("Cannot include more alpha spin-orbitals than there are")
        print ("Including the maximum amount")
        alpha_include = num_same_spin_orbitals - cbm_alpha_index - 1
        print ("New value for alpha_include = ", alpha_include)

    if cbm_beta_index + beta_include + 2 > 2*num_same_spin_orbitals:
        print ("Cannot include more beta spin-orbitals than there are")
        print ("Including the maximum amount")
        beta_include = 2*num_same_spin_orbitals - cbm_beta_index - 1
        print ("New value for beta_include = ", beta_include) 

    # Make ground state SD
    for i in range( cbm_alpha_index - alpha_include, cbm_alpha_index + 1):
        alpha_electrons.append( i )
    for i in range( cbm_beta_index  - beta_include, cbm_beta_index  + 1):
        beta_electrons.append( -i )

    gs_SD = alpha_electrons[:] + beta_electrons[:] 

    # Make excited state SDs
    es_SD    = []
    SD_basis = []
   
    # Excite only alpha electrons 
    if excitation_type == 0:
        for i in alpha_electrons:
            for j in range(cbm_alpha_index + 1, cbm_alpha_index + 2 + alpha_include):

                #es_sd = [ j if electron == i else electron for electron in gs_SD[:] ]  # compact version
                es_sd = []
                for electron in gs_SD[:]:
                    if electron == i:
                        es_sd.append( j ) 
                    else:
                        es_sd.append( electron )
                es_SD.append( es_sd )

    # Excite only beta electrons
    elif excitation_type == 1:   
        for i in beta_electrons:
            for j in range(cbm_beta_index + 1, cbm_beta_index + 2 + beta_include):

                #es_sd = [ -j if electron == i else electron for electron in gs_SD[:] ]  # compact version
                es_sd = []
                for electron in gs_SD[:]:
                    if electron == i:
                        es_sd.append( -j )
                    else:
                        es_sd.append( electron )
                es_SD.append( es_sd )
    
    # Excite both alpha and beta electrons
    else:
        for i in alpha_electrons:
            for j in range(cbm_alpha_index + 1, cbm_alpha_index + 2 + alpha_include):

                #es_sd = [ j if electron == i else electron for electron in gs_SD[:] ]  # compact version
                es_sd = []
                for electron in gs_SD[:]:
                    if electron == i:
                        es_sd.append( j )
                    else:
                        es_sd.append( electron )
                es_SD.append( es_sd )

        for i in beta_electrons:
            for j in range(cbm_beta_index + 1, cbm_beta_index + 2 + beta_include):

                #es_sd = [ -j if electron == i else electron for electron in gs_SD[:] ]  # compact version
                es_sd = []
                for electron in gs_SD[:]:
                    if electron == i:
                        es_sd.append( -j )
                    else:
                        es_sd.append( electron )
                es_SD.append( es_sd )

    SD_basis.append( gs_SD )
    for i in range(len(es_SD)):
        SD_basis.append( es_SD[i] )

    return SD_basis




def get_Lowdin(S):
    """  
    Find the S_i_half for the S matrix - alpha and beta components

    Args: 
        S ( CMATRIX(2N, 2N) ): is a matrix of MO overlaps. It has a block structure as:

            .. math::
                S = 
                \\begin{vmatrix}
                S_{aa}  & S_{ab} \\\
                S_{ba}  & S_{bb}
                \\end{vmatrix}

            Here, S_xy are the overlaps of the MOs for spin channels x and y (alpha, beta) - only
            spatial components of the orbitals are taken into account here.
            Here, N - is the total number of orbitals (double occupancies)        

    Returns:
        tuple: (S_aa_i_half, S_bb_i_half), where:

            * S_aa_i_half ( CMATRIX(N,N) ): S_aa^{-1/2} - inverse square root matrix for the alpha-alpha block
            * S_bb_i_half ( CMATRIX(N,N) ): S_bb^{-1/2} - inverse square root matrix for the beta-beta block
          
    """

    nstates = int(S.num_of_cols/2)  # division by 2 because it is a super-matrix

    alp = list(range(0,nstates))
    bet = list(range(nstates, 2*nstates))

    S_aa = CMATRIX(nstates, nstates)
    S_bb = CMATRIX(nstates, nstates)

    pop_submatrix(S, S_aa, alp, alp)
    pop_submatrix(S, S_bb, bet, bet)

    is_inv = FullPivLU_rank_invertible(S_aa)
    if is_inv[1] != 1:
        print("Error, S_aa is not invertible, Exiting Program");  sys.exit(0)

    is_inv = FullPivLU_rank_invertible(S_bb)
    if is_inv[1] != 1:
        print("Error, S_bb is not invertible, Exiting Program");  sys.exit(0)

    S_aa_half = CMATRIX(nstates,nstates)
    S_aa_i_half = CMATRIX(nstates,nstates)
    sqrt_matrix(S_aa, S_aa_half, S_aa_i_half)

    S_bb_half = CMATRIX(nstates,nstates)
    S_bb_i_half = CMATRIX(nstates,nstates)
    sqrt_matrix(S_bb, S_bb_half, S_bb_i_half)

    return S_aa_i_half, S_bb_i_half




def apply_normalization(S, St):
    """

    Transforms the input transition density matrix computed with potentially
    non-orthogonalized orbitals such that it would correspond to the properly
    orthonormalized ones

    Args: 
        S ( CMATRIX(2N, 2N) ): is a matrix of MO overlaps S_ij = <i|j>. It has a block structure as:

            .. math::
                S = 
                \\begin{vmatrix}
                S_{aa}  & S_{ab} \\\
                S_{ba}  & S_{bb}
                \\end{vmatrix}

            Here, S_xy are the overlaps of the MOs for spin channels x and y (alpha, beta) - only
            spatial components of the orbitals are taken into account here.
            Here, N - is the total number of orbitals (double occupancies)        

        St ( CMATRIX(2N, 2N) ): the transition density matrix St_ij = <i|d/dt|j>. It has a block structure as:

            .. math::
                St = 
                \\begin{vmatrix}
                St_{aa}  & St_{ab} \\\
                St_{ba}  & St_{bb}
                \\end{vmatrix}

            Here, St_xy are the transition density matrix for spin channels x and y (alpha, beta) - only
            spatial components of the orbitals are taken into account here.
            Here, N - is the total number of orbitals (double occupancies)        

    Returns:
        None: but the input matrix ```St``` is changed

    
    """

    nsteps = len(S)

    nstates = int(St[0].num_of_cols/2)  # division by 2 because it is a super-matrix
    
    alp = list(range(0,nstates))
    bet = list(range(nstates, 2*nstates))

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




def get_Lowdin_general(S):
    """  
    Find the S_i_half for the S matrix
    Args: 
        S ( CMATRIX(N, N) ): is a matrix of MO overlaps. 
    Returns:
        tuple: S_i_half, where:
            * S_i_half ( CMATRIX(N,N) ): S^{-1/2} - inverse square root matrix
          
    """

    nstates = int(S.num_of_cols)  # division by 2 because it is a super-matrix
    is_inv = FullPivLU_rank_invertible(S)
    if is_inv[1] != 1:
        print("Error, S is not invertible, Exiting Program");  sys.exit(0)
    S_half   = CMATRIX(nstates,nstates)
    S_i_half = CMATRIX(nstates,nstates)
    sqrt_matrix(S, S_half, S_i_half)
    return S_i_half




def apply_normalization_general(S, St):
    """

    Transforms the input transition density matrix computed with potentially
    non-orthogonalized orbitals such that it would correspond to the properly
    orthonormalized ones

    Args: 
        S  ( CMATRIX(N, N) ): is a matrix of MO overlaps S_ij = <i|j>
        St ( CMATRIX(N, N) ): the transition density matrix St_ij = <i|d/dt|j>

    Returns:
        None: but the input matricies ```S``` and ```St``` are changed   
    """

    nsteps  = len(S)
    nstates = int(St[0].num_of_cols)  # division by 2 because it is a super-matrix

    # For St
    for i in range(0, nsteps-1):
        U1 = get_Lowdin_general(S[i])   # time n
        U2 = get_Lowdin_general(S[i+1]) # time n+1          
        #print("St matrix before")
        #St[i].show_matrix()
        St_normalized = U1.H() * St[i] * U2
        #print("St matrix after")
        #St_normalized.show_matrix()
        push_submatrix(St[i], St_normalized, list(range(0,nstates)), list(range(0,nstates)))

    # For S
    for i in range(0, nsteps):
        U1 = get_Lowdin_general(S[i]) # time n
        #print("S matrix before")
        #S[i].show_matrix()
        S_normalized = U1.H() * S[i] * U1
        #print("S matrix after")
        #S_normalized.show_matrix()
        push_submatrix(S[i], S_normalized, list(range(0,nstates)), list(range(0,nstates)))





def make_cost_mat(orb_mat_inp, en_mat_inp, alpha):
    """

    Makes the cost matrix from a given TDM and information on states' energies

    Args:    
        orb_mat_inp  ( CMATRIX(nstates,nstates) or MATRIX(nstates,nstate) ): the transition density matrix
            TDM in a given basis. Here, ```nstates``` - the number of states (e.g. the number of doubly-occupied
            orbitals )

        en_mat_inp ( MATRIX(nstates, nstates) ): Matrix of energies in a given basis [units: a.u.]

        alpha ( float ): Parameter controlling the range of the orbitals that can participate in
            the reordering. Setting is to 0 makes all orbitals be considered for reordering
            Setting it to a large number makes the effective number of orbitals participating
            in the reordering smaller - this can be used to turn off the reordering. [units: a.u.^-1]

    Returns: 
        MATRIX(nstates, nstates): the matrix of the cost values for different pairs of states

    """

    nstates = orb_mat_inp.num_of_cols
    cost_mat = MATRIX(nstates, nstates)

    for a in range(0,nstates):
        for b in range(0,nstates):

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

    Args:
        St ( list of CMATRIX(nstates, nstates) ): TDM for each timestep
        E ( list of CMATRIX(nstates, nstates) ): energies of all states at every step
        params ( dictionary ): parameters controlling the reordering
            * **params["do_state_reordering"]** ( int ): option to select the state reordering algorithm 
                Available options:
                    - 1: older version developed by Kosuke Sato, may not the working all the times
                    - 2: Munkres-Kuhn (Hungarian) method [default]

            * **params["state_reordering_alpha"]** ( double ): a parameter that controls how 
                many states will be included in the reordering

    Returns:
        None: but changes the input St object

    """

    critical_params = [ ]
    default_params = { "do_state_reordering":2, "state_reordering_alpha":0.0 }
    comn.check_input(params, default_params, critical_params)

    nsteps = len(St)
    nstates = int(St[0].num_of_cols/2) # division by 2 because it is a super-matrix

    alp = list(range(0,nstates))
    bet = list(range(nstates, 2*nstates))

    # Initialize the cumulative permutation as the identity permutation
    perm_cum_aa = intList() # cumulative permutation for alpha spatial orbitals
    perm_cum_bb = intList() # cumulative permutation for beta  spatial orbtials
    for a in range(0,nstates):
        perm_cum_aa.append(a)
        perm_cum_bb.append(a)

    # Current permutation
    perm_t_aa = intList() 
    perm_t_bb = intList()
    for a in range(0,nstates):
        perm_t_aa.append(a)
        perm_t_bb.append(a)


    # Temporary matrices for the Hungarian method
    aa = CMATRIX(nstates, nstates); ab = CMATRIX(nstates, nstates)
    ba = CMATRIX(nstates, nstates); bb = CMATRIX(nstates, nstates)

    en_mat_aa = CMATRIX(nstates, nstates); 
    en_mat_bb = CMATRIX(nstates, nstates)


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




def do_phase_corr(cum_phase1, St, cum_phase2, phase_i):
    """

    This function changes the St matrix according to
    the previous cumulative phases and the current 
    phase correction as:

    St = <bra|ket>

    St -> St = F_n * St * (F_{n+1})^+ = F_n * St * (F_{n})^+ * (f_{n+1})^+

    Args:
        cum_phase1 ( CMATRIX(nstates, 1) ): cumulative phase corrections up to step n (F_n) for bra-vectors
        St         ( CMATRIX(nstates, nstates) ): input/output TDM to be processed: 
            could be alpha-alpha, beta-beta, alpha-beta, or beta-alpha sub-blocks
        cum_phase2 ( CMATRIX(nstates, 1) ): cumulative phase corrections up to step n (F_n) for ket-vectors
        phase_i    ( CMATRIX(nstates, 1) ): the current step phase corrections (f_{n+1}) for a given pair of vectors


    Returns: 
        None: but changes the input matrix St
  
    """
   
    nstates = St.num_of_rows

    ### Correct the TDM matrix ###
    for a in range(0,nstates):
        for b in range(0,nstates):
            fab = cum_phase1.get(a) * cum_phase2.get(b).conjugate() * phase_i.get(b).conjugate()
            St.scale(a,b, fab)



def apply_phase_correction(St):
    """Performs the phase correction according to:         
    Akimov, A. V. J. Phys. Chem. Lett, 2018, 9, 6096

    Args:
        St ( list of CMATRIX(N,N) ): St_ij[n] = <i(n)|j(n+1)> transition density matrix for 
            the timestep n, where N is the number of spin-orbitals in the active space. 
            Spin-orbitals, not just orbitals! So it is composed as:

            .. math::
                St = 
                \\begin{vmatrix}
                St_{aa}  & St_{ab} \\\
                St_{ba}  & St_{bb}
                \\end{vmatrix}

    Returns: 
        None: but changes the input St matrices

    """

    nsteps = len(St)
    nstates = int(St[0].num_of_cols/2)  # division by 2 because it is a super-matrix

    alp = list(range(0,nstates))
    bet = list(range(nstates, 2*nstates))

    ### Initiate the cumulative phase correction factors ###    
    cum_phase_aa = CMATRIX(nstates,1)  # F(n-1)  cumulative phase
    cum_phase_bb = CMATRIX(nstates,1)  # F(n-1)  cumulative phase
    for a in range(0,nstates):
        cum_phase_aa.set(a, 0, 1.0+0.0j)
        cum_phase_bb.set(a, 0, 1.0+0.0j)


    St_aa = CMATRIX(nstates, nstates); St_ab = CMATRIX(nstates, nstates)
    St_ba = CMATRIX(nstates, nstates); St_bb = CMATRIX(nstates, nstates)

    for i in range(0, nsteps):

        pop_submatrix(St[i], St_aa, alp, alp)
        pop_submatrix(St[i], St_bb, bet, bet)
        pop_submatrix(St[i], St_ab, alp, bet)
        pop_submatrix(St[i], St_ba, bet, alp)

        ### Compute the instantaneous phase correction factors for diag. blocks ###
        phase_i_aa = compute_phase_corrections(St_aa)   # f(i)
        phase_i_bb = compute_phase_corrections(St_bb)   # f(i)       

        ### Do the  phase correstions for the diag. blocks ###
        do_phase_corr(cum_phase_aa, St_aa, cum_phase_aa, phase_i_aa)
        do_phase_corr(cum_phase_bb, St_bb, cum_phase_bb, phase_i_bb)
        
        ### Do the  phase correstions for the off-diag. blocks ###
        do_phase_corr(cum_phase_aa, St_ab, cum_phase_bb, phase_i_bb)
        do_phase_corr(cum_phase_bb, St_ba, cum_phase_aa, phase_i_aa)

        ### Push the corrected diag. blocks to orig. St matrix ###
        push_submatrix(St[i], St_aa, alp, alp);   push_submatrix(St[i], St_ab, alp, bet)
        push_submatrix(St[i], St_ba, bet, alp);   push_submatrix(St[i], St_bb, bet, bet)

        ### Update the cumulative phase correction factors for diag. blocks ###
        for a in range(0,nstates):
            cum_phase_aa.scale(a, 0, phase_i_aa.get(a))
            cum_phase_bb.scale(a, 0, phase_i_bb.get(a))




def sac_matrices(coeff, basis, S_ks):
    """
    This function makes the Phi-to-Chi (P2C) transformation matrix.
    Normalization factros for the Chi states are computed based on the 
    overlaps of Phi states.
    < Chi_i | Chi_j > = 1
    = N_i * N_j * < Phi_i - Phi_i' | Phi_j - Phi_j' > = 1

    coeff [List of lists] - P2C as initialized by the user
    basis [Phi basis] - as initialized by the user 
    S_ks  [CMATRIX] - Time overlap matrix of elementary KS orbtials, from step2
    """

    n_chi = len(coeff)
    n_phi = len(coeff[0])

    P2C = CMATRIX(n_phi, n_chi)
    for j in range(0,n_chi):
        for i in range(0,n_phi):
            P2C.set(i,j,coeff[j][i]*(1.0+0.0j) )

    # Compute the overlaps of the SDs:
    #
    Ssd = mapping.ovlp_mat_arb(basis, basis, S_ks)

    # Normalize the Chi wavefunctions #
    norm = (P2C.H() * Ssd * P2C).real()
    for i in range(0,n_chi):
        if norm.get(i,i) > 0.0:
            P2C.scale(-1, i, 1.0/math.sqrt(norm.get(i,i)) )
        else:
            print("Error in CHI normalizaiton: some combination gives zero norm\n")
            sys.exit(0)

    return P2C


def scale_H_vib(hvib, en_gap, dNAC, sc_nac_method):
    """
    This function scales the energies and NACs in the vibrionic Hamiltonian
    in the Chi basis. 

    Args:   
        hvib ( list of CMATRIX objects ): 
            CMATRIXlist of vibronic hamiltonians in the Chi basis
        en_gap ( float ): The desired energy gap (E_1 - E_0), for the Chi basis
        dNAC ( list of lists of (list, float) ):
            The scaling terms by which specific nacs will 
            be scaled datatype = list of lists of (list, float)

            [  [ [i,j], val ], ...  ]          

            n and n+1 are the col (and thereby row) indicies of 
            the nacs to be scaled by the value val 
        sc_nac_method ( int ): The method used to scale NACs in the Chi basis, 
            chosen by the user.
            If sc_nac_method = 1, then the NACs are scaled by the ivnerse of the
            magnitude of the change in energy, according to Lin et al.

            Reference: Lin, Y. & Akimov, A. V. J. Phys. Chem. A (2016) 
    """

    traj_len = len(hvib)

    # Do the scaling
    nrows = hvib[0].num_of_rows
    ncols = hvib[0].num_of_cols   
    for i in range(0,traj_len):

        prev_gap = hvib[i].get(1,1) - hvib[i].get(0,0)
        shift = en_gap - prev_gap 

        # Scale the energy gap, by adding the scaling factors to all excited states
        # This is to keep the ground state enegy equal to 0.0 by definition
        for j in range(1,ncols):
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
            for j in range(1,ncols):
                scl_nac = prev_gap / (shift + prev_gap)
                hvib[i].scale(0, j, scl_nac)
                hvib[i].scale(j, 0, scl_nac)

    return hvib                    


def compute_Hvib(basis, St_ks, E_ks, dE, dt):
    """Compute the vibronic Hamiltonian matrix
    
    Args:    
        basis ( list of lists of integers ): defines the basis of Slater Determinants, 
            such that: basis[iSD][iks] is the indicator of the spin-orbital occupied by 
            the electron iks in the Slater Determinant iSD

            Example: 

                The following example defines a ground state SD (the lowest KS of the active space) and two 
                single excitations, which are different from each other by two spin flips of the electrons
                The convention is to start indexing from 1 (corresponds to index 0 in the KS matrices)
                Positive - for alpha electrons, negative - for beta electrons
                Need to be consistent: [ -1, 2 ] and [ 2, -1 ] are treated differently, this is needed for spin-adaptation

                >> basis = [ [ 1,-1 ], [ 1,-2 ], [ 2,-1 ] ]

                The next example is for a system of 4 electrons and hole excitations
                >> basis = [ [ 1,-1, 2, -2 ], [ 3, -1, 2, -2 ], [ 1, -3, 2, -2 ] ]

                                                                       
        St_ks ( CMATRIX(2*norbs, 2*norbs) ): transition density matrix in the KS spin-orbitals basis, where 
            norb - the number of double-occupied orbitals.

        E_ks ( CMATRIX(2*norbs, 2*norbs) ): the orbital energies in the KS spin-orbitals basis, where 
            norb - the number of double-occupied orbitals.

        dE ( list of doubles ): define corrections of the SD state energies in comparison to 
            the energy give by the sum energies of the occupied spin-orbitals.
            The convention is: dE[iSD] is the correction to energy of the SD with index iSD. 
            This is a constant correction - same for all energies in the set [units: Ha] 

            Example:
                For instance, for the SD examples above, the corrections could be:
                >> dE = [0.0, 0.01, 0.05]

        dt ( double ): the timestep for MD integrations [units: a.u.]

    Returns: 
        CMATRIX(nstates, nstates) The Vibronic Hamiltonian

    """

    St    = mapping.ovlp_mat_arb(basis, basis, St_ks) 
    H_el  = mapping.energy_mat_arb(basis, E_ks, dE)
    H_vib = H_el - (0.5j/dt)*CMATRIX((St-St.H()).real())

    return H_vib




def run(S_dia_ks, St_dia_ks, E_dia_ks, params):
    """
    The procedure to converts the results of QE calculations (KS orbital energies and
    time-overlaps = transition density matrices in the KS basis) to the generic Hvib matrices, 
    which (optionally) account for:   

    - enforces orthogonalization of the input KS states
    - state reordering
    - phase corrections
    - multi-electron wavefunction (Slater determinants) and spin-adaptation

    Args:
        S_dia_ks ( list of lists of CMATRIX objects ): overlaps of the KS orbitals along trajectories
            for each data set. Such that S_dia_ks[idata][istep].get(i,j) is <i(istep)|j(istep)> for the 
            trajectory (=dataset) ```idata```.

        St_dia_ks ( list of lists of CMATRIX objects ): time-overlaps (=transition density matrices) 
            in the basis of the KS orbitals along trajectories for each data set. 
            Such that St_dia_ks[idata][istep].get(i,j) is <i(istep)|j(istep+1)> for the trajectory (=dataset) ```idata```.

        E_dia_ks ( list of lists of CMATRIX objects ): energies the KS orbitals at the mid-points along trajectories
            for each data set. Such that E_dia_ks[idata][istep].get(i,i) is 0.5*(E_i(istep) + E_i(istep+1)) for the 
            trajectory (=dataset) ```idata```

        params ( dictionary ): Control paramerter of this type of simulation. Can include the follwing keys:

            * **params["SD_basis"]** ( list of lists of ints ): define the Slater Determinants basis
                The convention is:  params["SD_basis"][iSD][iks] is the indicator of the spin-orbital occupied by 
                the electron iks in the Slater Determinant iSD [required!]

                Example: 

                    The following example defines a ground state SD (the lowest KS of the active space) and two 
                    single excitations, which are different from each other by two spin flips of the electrons
                    The convention is to start indexing from 1 (corresponds to index 0 in the KS matrices)
                    Positive - for alpha electrons, negative - for beta electrons
                    Need to be consistent: [ -1, 2 ] and [ 2, -1 ] are treated differently, this is needed for spin-adaptation

                    >> params["SD_basis"] = [ [ 1,-1 ], [ 1,-2 ], [ 2,-1 ] ]

                    The next example is for a system of 4 electrons and hole excitations
                    >> params["SD_basis"] = [ [ 1,-1, 2, -2 ], [ 3, -1, 2, -2 ], [ 1, -3, 2, -2 ] ]


            * **params["SD_energy_corr"]** ( list of doubles ): define corrections of the SD state energies in comparison to 
                the energy give by the sum energies of the occupied spin-orbitals.
                The convention is: params["SD_energy_corr"][iSD] is the correction to energy of the SD with index iSD. 
                This is a constant correction - same for all energies in the set [units: Ha] [required!]

                Example:
                    For instance, for the SD examples above, the corrections could be:
                    >> params["SD"] = [0.0, 0.01, 0.05]
                                     
            * **params["CI_basis"]** ( list of lists of complex number ): configuration interaction coefficients 
                that define a superpositions to SDs that are considered the states of interest, e.g. spin-adapted configurations
                The convention is: params["CI_basis"][iCI][iSD] is a coefficient of ```iSD```-th SD in the expansion of the CI
                with index ```iCI```. These coefficients don't have to account for the overal CI's normalization - the 
                normalization will be done on the go. [required!]

                Example:

                    For the SD example above we can construct the following combinations:
                    >> params["CI_basis"] = [ [1.0, 0.0, 0.0 ], 
                                              [0.0, 1.0,-1.0 ],
                                              [0.0, 1.0, 1.0 ] 
                                            ]

            * **params["output_set_paths"]** ( list of strings ): the directory pathes where the resulting files 
                are to be written (if so!). If you don't plan on writing the files, just provide a list of empty strings
                or whatever else - they will not be used in that case. The number of the strings should be equal to 
                the number of the input data sets, e.g. to len(St_dia_ks)  [required!]

                  
            * **params["dt"]** ( double ): nuclear dynamics integration time step [units: a.u. of time, default: 41.0]
 
            * **params["do_orthogonalization"]** ( int ): the option to do Lowdin orthogonalization of the orbitals - using 
                the "raw" overlaps (at the same time). This option is needed because the wavefunction output by QE are 
                not exactly orthonormal (because of the use of pseudopotentials). So before we use them (implicitly) 
                in the rest of the calculations here, we may need to account for this non-ideality effect.
                Options:
        
                - 0: don't do the orthogonalization - this is the same as in Pyxaid [default]
                - 1: do the orthogonalization

            * **params["do_state_reordering"]** ( int ): the option to control the state reordering
                Options:

                - 0: no state reordering - same as in Pyxaid
                - 1: older method (is not robust, may or may not work) 
                - 2: Hungarian algorithm [default]

            * **params["state_reordering_alpha"]** ( double ): the parameter that controls the width of 
                the energy interval within wich the state reordering is in effect. Zero value means all 
                available orbitals, larger positive value decreases the width of the window. This parameter
                is not in effect unless the Hungarian algorithm is selected [default: 0.0]

            * **params["do_phase_correction"]** ( int ): option to do the phase correction

                - 0 - don't do 
                - 1 - do it [default]

            * **params["do_output"]** ( int ): whether to print out the Hvib matrices ( = 1) to the files or not ( = 0).

            * **params["Hvib_re_prefix"]** ( string ): common prefix of the output files with real part of the vibronic 
                Hamiltonian at all times [default: "Hvib_"]

            * **params["Hvib_re_suffix"]** ( string ): common suffix of the output files with real part of the vibronic 
                Hamiltonian at all times [default: "_re"]

            * **params["Hvib_im_prefix"]** ( string ): common prefix of the output files with imaginary part of the vibronic 
                Hamiltonian at all times [default: "Hvib_"]

            * **params["Hvib_im_suffix"]** ( string ): common suffix of the output files with imaginary part of the vibronic 
                Hamiltonian at all times [default: "_im"]

    Returns:
        list of lists of CMATRIX(N,N): Hvib, such that:
            Hvib[idata][istep] is a CMATRIX(N,N) containing the vibronic Hamiltonian for the 
            trajectory (dataset) ```idata``` and for the timestep ```istep```. Here, N is the number
            of states included in the active space.

    """


    #====== Defaults and local parameters ===============

    critical_params = [ "SD_basis", "SD_energy_corr", "CI_basis", "output_set_paths" ]
    default_params = { "dt":1.0*units.fs2au, 
                       "do_orthogonalization":0,
                       "do_state_reordering":2, "state_reordering_alpha":0.0,
                       "do_phase_correction":1,
                       "do_output":0,
                       "Hvib_re_prefix":"Hvib_", "Hvib_im_prefix":"Hvib_",
                       "Hvib_re_suffix":"_re", "Hvib_im_suffix":"_im",

                     }
    comn.check_input(params, default_params, critical_params)
 
    dt = params["dt"]
    do_orthogonalization = params["do_orthogonalization"]
    do_phase_correction = params["do_phase_correction"]
    do_state_reordering = params["do_state_reordering"]
    ndata =  len(St_dia_ks)
    nsteps = len(St_dia_ks[0])
    nstates = len(params["CI_basis"])  # the number of CI states to consider

    do_output = params["do_output"]

    #====== Generic sanity checks ===============
    if do_output:
        if(ndata) != len(params["output_set_paths"]):
            print("Error: The number of output sets paths should be the same as the number of data sets\n")
            print("ndata = ", ndata)
            print("len(params[\"output_set_paths\"]) = ", len(params["output_set_paths"]))
            print("Exiting...\n")
            sys.exit(0)



    #====== Calculations  ===============
    H_vib = []

    for idata in range(0,ndata):

        # 1. Do the KS orbitals orthogonalization 
        if do_orthogonalization > 0:
            apply_normalization(S_dia_ks[idata], St_dia_ks[idata])

        # 2. Apply state reordering to KS
        if do_state_reordering > 0:
            apply_state_reordering(St_dia_ks[idata], E_dia_ks[idata], params)

        # 3. Apply phase correction to KS
        if do_phase_correction > 0:
            apply_phase_correction(St_dia_ks[idata])
       
        
        Hvib = []
        for i in range(0,nsteps):
        
            # 4. Construct the Hvib in the basis of Slater determinants (SDs)
            hvib_sd = compute_Hvib(params["SD_basis"], St_dia_ks[idata][i], E_dia_ks[idata][i], params["SD_energy_corr"], dt) 
      
            # 5. Convert the Hvib to the basis of symmery-adapted configurations (SAC)
            SD2CI = sac_matrices(params["CI_basis"], params["SD_basis"], S_dia_ks[idata][i])
            hvib_ci = SD2CI.H() * hvib_sd * SD2CI
            Hvib.append( hvib_ci )


        if do_output:
            # Output the resulting Hamiltonians
            for i in range(0,nsteps):
                re_filename = params["output_set_paths"][idata] + params["Hvib_re_prefix"] + str(i) + params["Hvib_re_suffix"]
                im_filename = params["output_set_paths"][idata] + params["Hvib_im_prefix"] + str(i) + params["Hvib_im_suffix"]        
                Hvib[i].real().show_matrix(re_filename)
                Hvib[i].imag().show_matrix(im_filename)


        H_vib.append(Hvib)        
        
    return H_vib




def map_Hvib(H, basis, dE):

    nbasis = -1 # doesn't matter

    H_vib  = mapping.energy_mat_arb(basis, H, dE)

    nstates = len(basis)
    for i in range(0,nstates):
        for j in range(0,nstates):

            res = delta(Py2Cpp_int(basis[i]), Py2Cpp_int(basis[j]) )

            if res[0] != 0:
                if res[1] * res[2] > 0:
                    a = mapping.sd2indx([res[1], res[2]], nbasis, False) 
                    H_vib.add(i, j, H.get(a[0], a[1]) )

    return H_vib



def pyxaid2libra(Hvib_pyxaid, params):
    """
    """

    #====== Defaults and local parameters ===============

    critical_params = [ "SD_basis", "SD_energy_corr", "CI_basis", "output_set_paths" ]
    default_params = { "do_output":0,
                       "Hvib_re_prefix":"Hvib_", "Hvib_im_prefix":"Hvib_",
                       "Hvib_re_suffix":"_re", "Hvib_im_suffix":"_im",
                     }
    comn.check_input(params, default_params, critical_params)

      
    ndata =  len(Hvib_pyxaid)
    nsteps = len(Hvib_pyxaid[0])
    norbs = Hvib_pyxaid[0][0].num_of_cols
    nstates = len(params["CI_basis"])  # the number of CI states to consider
    do_output = params["do_output"]


    S = CMATRIX(2*norbs, 2*norbs)
    H = CMATRIX(2*norbs, 2*norbs)
    s = CMATRIX(norbs, norbs)
    s.identity()
    alp = list(range(0, norbs))
    bet = list(range(norbs, 2*norbs))

    push_submatrix(S, s, alp, alp)
    push_submatrix(S, s, alp, bet)
    push_submatrix(S, s, bet, alp)
    push_submatrix(S, s, bet, bet)

    #====== Generic sanity checks ===============
    if do_output:
        if(ndata) != len(params["output_set_paths"]):
            print("Error: The number of output sets paths should be the same as the number of data sets\n")
            print("ndata = ", ndata)
            print("len(params[\"output_set_paths\"]) = ", len(params["output_set_paths"]) )
            print("Exiting...\n")
            sys.exit(0)



    #====== Calculations  ===============
    H_vib = []

    for idata in range(0,ndata):

        Hvib = []
        for i in range(0,nsteps):
        
            # 1. Construct the Hvib in the basis of Slater determinants (SDs)            
            push_submatrix(H, Hvib_pyxaid[idata][i], alp, alp)
            push_submatrix(H, Hvib_pyxaid[idata][i], alp, bet)
            push_submatrix(H, Hvib_pyxaid[idata][i], bet, alp)
            push_submatrix(H, Hvib_pyxaid[idata][i], bet, bet)
            hvib_sd = map_Hvib(H, params["SD_basis"], params["SD_energy_corr"]) 
      
            # 2. Convert the Hvib to the basis of symmery-adapted configurations (SAC)
            SD2CI = sac_matrices(params["CI_basis"], params["SD_basis"], S)
            hvib_ci = SD2CI.H() * hvib_sd * SD2CI
            Hvib.append( hvib_ci )


        if do_output:
            # Output the resulting Hamiltonians
            for i in range(0,nsteps):
                re_filename = params["output_set_paths"][idata] + params["Hvib_re_prefix"] + str(i) + params["Hvib_re_suffix"]
                im_filename = params["output_set_paths"][idata] + params["Hvib_im_prefix"] + str(i) + params["Hvib_im_suffix"]        
                Hvib[i].real().show_matrix(re_filename)
                Hvib[i].imag().show_matrix(im_filename)

        H_vib.append(Hvib)        
        
    return H_vib



def apply_state_reordering_general(St, E, params):
    """
    Performs the state's identity reordering in a given basis for all time steps.
    This is reflects in the corresponding changess of the TDM.

    This function is for dat NOT in spin-block format

    Args:
        St ( list of CMATRIX(nstates, nstates) ): TDM for each timestep
        E ( list of CMATRIX(nstates, nstates) ): energies of all states at every step
        params ( dictionary ): parameters controlling the reordering
            * **params["do_state_reordering"]** ( int ): option to select the state reordering algorithm 
                Available options:
                    - 1: older version developed by Kosuke Sato, may not the working all the times
                    - 2: Munkres-Kuhn (Hungarian) method [default]
            * **params["state_reordering_alpha"]** ( double ): a parameter that controls how 
                many states will be included in the reordering
    Returns:
        None: but changes the input St object
    """

    critical_params = [ ]
    default_params = { "do_state_reordering":2, "state_reordering_alpha":0.0 }
    comn.check_input(params, default_params, critical_params)

    nsteps  = len(St)
    nstates = St[0].num_of_cols 

    # Initialize the cumulative permutation as the identity permutation
    perm_cum = intList() # cumulative permutation for alpha spatial orbitals
    for i in range(0,nstates):
        perm_cum.append(i)

    # Current permutation
    perm_t = intList() 
    for i in range(0,nstates):
        perm_t.append(i)
 
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

            # Permute rows 
            St[i].permute_rows(perm_t)

            # compute the cost matrices
            cost_mat = make_cost_mat(St[i], E[i], params["state_reordering_alpha"])

            # Solve the optimal assignment problem for diagonal blocks
            res = hungarian.maximize(cost_mat)

            # Convert the list of lists into the permutation object
            for row in res:
                perm_t[row[0]] = row[1]  # for < i | i > this becomes a new value: perm_t = P_{n+1}

            # Permute the blocks by col
            St[i].permute_cols(perm_t)



def apply_phase_correction_general(St):
    """Performs the phase correction according to:         
    Akimov, A. V. J. Phys. Chem. Lett, 2018, 9, 6096

    This function is for dat NOT in spin-block format

    Args:
        St ( list of CMATRIX(N,N) ): St_ij[n] = <i(n)|j(n+1)> transition density matrix for 
            the timestep n, where N is the number of states in the active space. 
            Spin-orbitals, not just orbitals! So it is composed as:

    Returns: 
        None: but changes the input St matrices
    """

    nsteps  = len(St)
    nstates = St[0].num_of_cols

    ### Initiate the cumulative phase correction factors ###    
    cum_phase = CMATRIX(nstates,1)  # F(n-1)  cumulative phase
    for i in range(0,nstates):
        cum_phase.set(i, 0, 1.0+0.0j)

    print("number of steps , nsteps= ", nsteps)

    for i in range(0, nsteps):

        ### Compute the instantaneous phase correction factors for diag. blocks ###
        phase_i = compute_phase_corrections(St[i]) # f(i)

        phase_i.show_matrix()

        ### Do the  phase corrections ###
        do_phase_corr(cum_phase, St[i], cum_phase, phase_i)

        ### Update the cumulative phase correction factors for diag. blocks ###
        for j in range(0,nstates):
            cum_phase.scale(j, 0, phase_i.get(j))




