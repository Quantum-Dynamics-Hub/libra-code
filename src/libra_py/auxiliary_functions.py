#*********************************************************************************
#*
#* Copyright (C) 2020 Mohammad Shakiba, Brendan Smith, Alexey V. Akimov
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: auxiliary_functions
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing with cube files
.. moduleauthors:: 
       Mohammad Shakiba, Brendan Smith, Alexey V. Akimov 
  
"""



import numpy as np
import math
import sys
import os

import time
import copy

import multiprocessing as mp

from liblibra_core import *

from libra_py import data_conv
from libra_py import cube_file_methods
from libra_py import CP2K_methods
from libra_py.workflows.nbra import mapping
from libra_py.workflows.nbra import step3
from libra_py import units

# Concurrency for large systems
import concurrent.futures








# This file is temp only. These functions will eventually be placed in Libra somewhere ...
def curr_and_final_step_job( istep, fstep, njobs, njob ):
    """
    This function is used to determine the initial and final step of a job when distributing
    the molecular dynamics trajectory over a number of jobs.
    
    Args:
    
        istep (int): The initial time step for molecular dynamics trajectory which 
	             the user wants to start the calculations with. This paramter 
		     is set in the input file for submitting the jobs.
        
	fstep (int): The final time step for molecular dynamics trajectory which 
	             the user wants to start the calculations with. This paramter 
		     is set in the input file for submitting the jobs.
		     
        njobs (int): The number of jobs specified by user.
	
	njob (int): The job number.
	
    Returns:
    
        job_init_step (int): The job initial time step.
	
	job_final_step (int): The job final time step.
	
    """

    total_steps = fstep - istep + 1

    nsteps_per_job = int(total_steps/njobs)

    if njob == 0:
        job_init_step  = istep
        job_final_step = istep + nsteps_per_job - 1

    elif njob > 0 and njob < ( njobs - 1 ):
        job_init_step  = istep + njob * nsteps_per_job - 1
        job_final_step = istep + ( njob + 1 ) * nsteps_per_job - 1

    elif njob == ( njobs - 1 ):
        job_init_step  = istep + njob * nsteps_per_job - 1
        job_final_step = fstep

    return job_init_step, job_final_step





def normalize_ci_coefficients(ci_coefficients_raw_unnorm):
    """
    This funciton normalizes the list of configuration interaction (CI) coefficients.
    
    Args:
    
        ci_coefficients_raw_unnorm (list): The list containing the lists of unnormalized CI coefficients.
	
    Returns:
    
        ci_coefficients_raw_norm (list): The list containing the lists of normalized CI coefficients.
    
    """
    # Number of states contributing in the excited state
    nstates = len(ci_coefficients_raw_unnorm)
    # Creating an empty list to store the normalized CI coefficients
    ci_coefficients_raw_norm = []
       
    for i in range(nstates):
        
	#### ordinary way without using numpy
        # Set up an initial parameter to compute the norm of the vector
	norm = 0.0
	# For each list of CI coefficients
        ci_coefficients_raw_norm.append( [] )
        
        for j in ci_coefficients_raw_unnorm[i]:
            # sum of CI coefficients square
            norm += j*j
        # Compute the norm by taking the square root of the sum to compute the norm
        norm = math.sqrt(norm)
        
        # numpy way to compute the normalized CI coefficients
        ci_coefficients_raw_unnorm[i] = np.array(ci_coefficients_raw_unnorm[i])
        ci_coefficients_raw_norm[i] = abs ( ci_coefficients_raw_unnorm[i] / norm )
        ci_coefficients_raw_norm[i] = list(ci_coefficients_raw_norm[i])
        
    
    return ci_coefficients_raw_norm





def get_es_output(params):
    """
    This function reads the information of the excited states from the log file of the single point calculations.
    
    Args:
    
        params (dict):
	
	    logfile_directory (str): The log files directory.
	    
	    es_software (str): The name of the software used to calculate the energy calculations.
	    
	    curr_step (int): The current time step of the calculations.
	    
    Returns:
    
        excitation_energies (numpy array): The excitation energies of the curr_step.
	
	ci_basis_raw (numpy array): The excited states which contains the occupied and virtual orbitals.
	
	ci_coefficients_raw_unnorm (numpy array): The excited states CI coefficients.
	
	spin_components (numpy array): Contains the excited states spin components ('alp' for alpha spin and 'bet' for beta spin)
	
    """
    
    logfile_directory = params["logfile_directory"]
    es_software = params["es_software"]
    curr_step = params["curr_step"]

    if es_software == "cp2k":
        
        logfile_name = logfile_directory+'/step_'+str(curr_step)+'.log'
        params.update({"logfile_name":logfile_name}) 
        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = \
        CP2K_methods.read_cp2k_tddfpt_log_file( params )
  
    return excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components





def integrate_cube_set( cubefiles_curr, cubefiles_prev, dv ):
    """
    This function 
    """

    overlap_matrix = np.zeros( ( len( cubefiles_prev ) , len( cubefiles_curr ) ) )

    for i in range( 0, len( cubefiles_prev ) ):
        for j in range( 0, len( cubefiles_curr ) ):
            overlap_matrix[i][j] = cube_file_methods.integrate_cube( cubefiles_prev[i], cubefiles_curr[j], dv )

    return overlap_matrix





def compute_cube_ks_overlaps( cubefiles_prev, params):
    """
    This function computes
    """

    curr_step = int(params["curr_step"])
    isUKS  = int(params["isUKS"])
    nprocs = int(params["nprocs"])
    
    print('---------------------------------------------------------')
    print('Starting the calculations by reading the cube files    \n')
    print('                    Step time %i\n'%(curr_step))
    print('-------------------------------------------------------\n')
    print('Reading cube files with %i number of processors'%nprocs)
    # First read all the .cube files (This is the most time consuming part)
    # Creating a pool and assigning the number of processors for mp.Pool
    pool = mp.Pool(processes=nprocs)
    
    # for time t
    print("\nWe are reading the cubes for step ", curr_step)
    cubefile_names_curr = CP2K_methods.cube_file_names_cp2k( params )
    print("\nIn compute_cube_ks_overlaps, the cubfile_names_curr are")
    print("\n", cubefile_names_curr)

    cubefiles_curr = []    
    cubefiles_curr = pool.map( cube_file_methods.read_cube, cubefile_names_curr ) 
    pool.close()
    
    # Calculate the dv element for integration
    dv = cube_file_methods.grid_volume( cubefile_names_curr[0] )

    # Compute the overlaps
    if isUKS == 1:

        # The cube files for alpha and beta spin for time t-1
        alp_cubes_prev = cubefiles_prev[0::2]
        bet_cubes_prev = cubefiles_prev[1::2]

        # The cube files for alpha and beta spin for time t
        alp_cubes_curr = cubefiles_curr[0::2]
        bet_cubes_curr = cubefiles_curr[1::2]

        zero_mat_alp = np.zeros( ( len( alp_cubes_prev ), len( alp_cubes_curr ) ) )
        zero_mat_bet = np.zeros( ( len( bet_cubes_prev ), len( bet_cubes_curr ) ) )

        """
        print('---------------------------------')
        print('Calculating the time overlap matrix...')
        # The time overlap between cube files of time t-1 and t
        St_alp_alp = integrate_cube_set( alp_cubes_prev, alp_cubes_curr, dv ) 
        St_bet_bet = integrate_cube_set( bet_cubes_prev, bet_cubes_curr, dv )
	
        # The overlap between cube files at times t-1 and t
        S_alp_alp_prev = integrate_cube_set( alp_cubes_prev, alp_cubes_prev, dv )
        S_bet_bet_prev = integrate_cube_set( bet_cubes_prev, bet_cubes_prev, dv )
        S_alp_alp_curr = integrate_cube_set( alp_cubes_curr, alp_cubes_curr, dv )
        S_bet_bet_curr = integrate_cube_set( bet_cubes_curr, bet_cubes_curr, dv )
	"""

        ### For using concurrency 
        with concurrent.futures.ThreadPoolExecutor(max_workers=nprocs) as executor:
            int_1 = executor.submit(integrate_cube_set, alp_cubes_prev, alp_cubes_curr, dv )
            int_2 = executor.submit(integrate_cube_set, bet_cubes_prev, bet_cubes_curr, dv )
            int_3 = executor.submit(integrate_cube_set, alp_cubes_prev, alp_cubes_prev, dv )
            int_4 = executor.submit(integrate_cube_set, bet_cubes_prev, bet_cubes_prev, dv )
            int_5 = executor.submit(integrate_cube_set, alp_cubes_curr, alp_cubes_curr, dv )
            int_6 = executor.submit(integrate_cube_set, bet_cubes_curr, bet_cubes_curr, dv )
        St_alp_alp = int_1.result()
        St_bet_bet = int_2.result()
        ### End using concurrency 
	
        # The overlap between cube files at times t-1 and t
        S_alp_alp_prev = int_3.result()
        S_bet_bet_prev = int_4.result()
        S_alp_alp_curr = int_5.result()
        S_bet_bet_curr = int_6.result()

    else:
        """
        # Spin non-polarized case
        St = np.zeros( ( len( cubefiles_prev ), len( cubefiles_curr ) ) )
        S_prev = np.zeros( ( len( cubefiles_prev ), len( cubefiles_prev ) ) )
        S_curr = np.zeros( ( len( cubefiles_curr ), len( cubefiles_curr ) ) )
    
        print('---------------------------------')
        print('Calculating the time overlap matrix...')
        for i in range( 0, len( cubefiles_prev ) ):
            for j in range( 0, len( cubefiles_curr ) ):
                S_prev[i][j] = cube_file_methods.integrate_cube( cubefiles_prev[i][0], cubefiles_prev[j][0], dv )
                S_curr[i][j] = cube_file_methods.integrate_cube( cubefiles_curr[i][0], cubefiles_curr[j][0], dv )
                St[i][j]     = cube_file_methods.integrate_cube( cubefiles_prev[i][0], cubefiles_curr[j][0], dv )
       """

        ### Using the concurrency
        with concurrent.futures.ThreadPoolExecutor(max_workers=nprocs) as executor:
            int_1 = executor.submit(integrate_cube_set, cubefiles_prev, cubefiles_prev, dv )
            int_2 = executor.submit(integrate_cube_set, cubefiles_curr, cubefiles_curr, dv )
            int_3 = executor.submit(integrate_cube_set, cubefiles_prev, cubefiles_curr, dv )
        S_prev = int_1.result()
        S_curr = int_2.result()
        St     = int_3.result()
        ### End using concurrency

        zero_mat_alp = np.zeros( ( len( S_prev ), len( S_curr ) ) )
        zero_mat_bet = np.zeros( ( len( S_prev ), len( S_curr ) ) )
        S_alp_alp_prev, S_bet_bet_prev = S_prev, S_prev
        S_alp_alp_curr, S_bet_bet_curr = S_curr, S_curr
        St_alp_alp, St_bet_bet = St, St
    
    S_ks_prev = data_conv.form_block_matrix( S_alp_alp_prev, zero_mat_alp, zero_mat_bet, S_bet_bet_prev )
    S_ks_curr = data_conv.form_block_matrix( S_alp_alp_curr, zero_mat_alp, zero_mat_bet, S_bet_bet_curr )
    St_ks     = data_conv.form_block_matrix( St_alp_alp,     zero_mat_alp, zero_mat_bet, St_bet_bet )

    return cubefiles_curr, S_ks_prev, S_ks_curr, St_ks 






def reindex_cpk2_sd_states( ks_orbital_homo_index, ks_orbital_indicies, sd_basis_states ):
    """

    sd_basis_states( list of lists of lists ): A list of Slater determinants, where each slater determinant is a excitation in the Kohn-Sham
                                               basis. This function assumes that all Kohn-Sham excitations are for alpha electrons. To
                                               differentiate between alpha and beta excitations, elements of sd_basis_states contain spin
                                               information.
               
                                               Ex) sd_basis_states[0] = [ [9,10], "alp" ] 
                                                   sd_basis_states[1] = [ [9,10], "bet" ]

    Returns:

         Reindexed Slater determinant basis states in terms of the ks_orbital_homo_index. 

         Ex.)   ks_orbital_homo_index = 9,    ks_orbital_indicies = [5,6,7,8,9,10,11,12,13,14]
                [ [9,10], "alp" ] --> [5,6]
                [ [9,10], "bet" ] --> [-15,-16]

         In general form, the reindexing take the form of:

                [ [9,10], "alp" ] --> [ index from 1 of 9 in ks_orbital_indicies, index from 1 of 10 in ks_orbital_indicies  ]
                [ [9,10], "bet" ] --> [ - ( index from 1 of 9 in ks_orbital_indicies + len(ks_orbital_indicies) ), - ( index from 1 of 10 in ks_orbital_indicies + len(ks_orbital_indicies) ) ]

    """

    # We need to update the indexing of the sd_basis - in terms of the rows and cols of St_KS
    # reindex ks orbs according to the matrix size
    n_alp_ks_orbs = len(ks_orbital_indicies)
    alp_homo_matrix_index = 0
    for i in range( n_alp_ks_orbs ):
        if ks_orbital_indicies[i] == ks_orbital_homo_index:
            alp_homo_matrix_index = i+1

    ks_orbs_new_index = [ i+1 for i in range(n_alp_ks_orbs) ]

    # Form excited state SDs
    excitations = []
    # For each Slater determinant basis state, which could have spin_component "alp" or "bet"
    for j in range( len( sd_basis_states ) ):

        if sd_basis_states[j][1] == "alp":
            initial_ks_orb = int( sd_basis_states[j][0][0] ) - ks_orbital_homo_index + alp_homo_matrix_index
            final_ks_orb   = int( sd_basis_states[j][0][1] ) - ks_orbital_homo_index + alp_homo_matrix_index
            excitations.append( [initial_ks_orb, final_ks_orb]  )
    
        elif sd_basis_states[j][1] == "bet":
            initial_ks_orb = int( sd_basis_states[j][0][0] ) - ks_orbital_homo_index + alp_homo_matrix_index + n_alp_ks_orbs
            final_ks_orb   = int( sd_basis_states[j][0][1] ) - ks_orbital_homo_index + alp_homo_matrix_index + n_alp_ks_orbs
            excitations.append( [-initial_ks_orb, -final_ks_orb]  )
    print( "excitations = ", excitations )
    #sys.exit(0)

    sd_basis = [ [] ]
    # We are now going to form the ground state Slater determinant. We first need a loop over all of the Kohn-Sham orbitals
    # beginning from index 1.
    # At this point we have a bunch of Slater determinants in the form of reindexed excitations in the variable "excitations" 
    for i in range( 1, 2*len( ks_orbital_indicies )-1 ):
        # Form ground state SD. Following the example in the function's documentation, the ground state Slater determinant should be:
        # [ 1, 2, 3, 4, 5, -11, -12, -13, -14,- 15 ]. Later on, we have the option of using only [ 5, -15 ] when computing the overlaps
        if i < alp_homo_matrix_index + 1:
            sd_basis[0].append( i )
        if i > n_alp_ks_orbs and i < n_alp_ks_orbs + alp_homo_matrix_index + 1:
            sd_basis[0].append( -i )
    #print ( sd_basis )

    # Now that we have done the ground state slater 
    for j in range( len( excitations ) ):
        sd_excitation = [ excitations[j][1] if x == excitations[j][0] else x for x in sd_basis[0] ]
        sd_basis.append( sd_excitation )
    print ( sd_basis )
    
    return sd_basis



def apply_state_reordering_ci(St, E, params):
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

            # Permute rows for this time-step
            St[i].permute_rows(perm_t)

            # compute the cost matrices
            cost_mat = make_cost_mat(St[i], E[i], params["state_reordering_alpha"])

            # Solve the optimal assignment problem for this time-step
            res = hungarian.maximize(cost_mat)

            # Convert the list of lists into the permutation object
            for row in res:
                perm_t[row[0]] = row[1]  # for < i | i > this becomes a new value: perm_t = P_{n+1}

            # Permute the blocks by col
            St[i].permute_cols(perm_t)


		
def form_Hvib_real( cp2k_log_file_name: str, time: int, min_band: int, max_band: int, isUKS: int):
    """
    This function forms the real part of the vibronic Hamiltonian by inserting the 
    energies on the diagonal of a zero matrix. The Hvib is in two-spinor format 
    containing the alpha and beta states energies.
    
    Args:
        
        cp2k_log_file_name (str): The output file produced by CP2K
        
        time (int): The time step of the molecular dynamics. For single point calculations
                    it should be set to 0.
                    
        min_band (int): The minimum state to be considered.
        
        max_band (int): The maximum state to be considered.
        
        isUKS (int): This flag is used whenever the UKS calculations were called in CP2K input. if it is set to
                     1 then UKS calculations were called and if not it will consider only alpha energies.
        
    Returns:
    
        Hvib_ks_re (2D numpy array): The diagonal matrix containing the energies of alpha and beta spins if
                                     UKS is set to True and the energies of alpha spin in a block form matrix.
        
    """

    # If the UKS calculations were set to True in CP2K input
    if isUKS == 1:
        # Read energies with alpha spin
        spin = 1
        E_alpha, total_energy = read_energies_from_cp2k_log_file( cp2k_log_file_name, time, min_band, max_band, spin )
        # Read energies with beta spin
        spin = 2
        E_beta, total_energy = read_energies_from_cp2k_log_file( cp2k_log_file_name, time, min_band, max_band, spin )
        # Now forming the diagonal matrix containing the Kohn-Sham energies of the alpha and beta spins
        Hvib_ks_re = np.diag( np.concatenate( ( E_alpha, E_beta ) ) )
    # If there is no UKS calculations set
    else:
        spin = 1
        E_alpha, total_energy = read_energies_from_cp2k_log_file( cp2k_log_file_name, time, min_band, max_band, spin )
        Hvib_ks_re = np.diag( np.concatenate( ( E_alpha, E_alpha ) ) )
    
    return Hvib_ks_re, total_energy



		
