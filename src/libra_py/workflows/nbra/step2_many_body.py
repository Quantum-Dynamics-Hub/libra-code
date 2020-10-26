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
.. module:: step2_many_body
   :platform: Unix, Windows
   :synopsis: This module implements functions for computing the overlap and nonadiabatic couplings matrices.
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

# Concurrency for large systems
import concurrent.futures

from liblibra_core import *

from libra_py import data_conv
from libra_py import cube_file_methods
from libra_py import CP2K_methods
from libra_py import Gaussian_methods
from libra_py.workflows.nbra import mapping
from libra_py.workflows.nbra import step3
from libra_py import units

import util.libutil as comn


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
        ci_coefficients_raw_norm[i]   = ci_coefficients_raw_unnorm[i] / norm
        ci_coefficients_raw_norm[i]   = list(ci_coefficients_raw_norm[i]) 
    
    return ci_coefficients_raw_norm





def get_excitation_analysis_output(params):
    """
    This function reads the information of the excited states from the log file of the single point calculations.
    
    Args:
    
        params (dict):

            logfile_directory (str): The log files directory.

            es_software (str): The name of the software used to calculate the energy calculations.
        
            curr_step (int): The current time step of the calculations.
            
            isUKS (int): This parameter is set for spin restricted and unrestricted calculations. When it is
                         set to 1 it means that unrestricted calculations were set in the input file otherwise 
                         it is restricted.
    
    Returns:
    
        excitation_energies (numpy array): The excitation energies of the curr_step.

        ci_basis_raw (numpy array): The excited states which contains the occupied and virtual orbitals.

        ci_coefficients_raw_unnorm (numpy array): The excited states CI coefficients.

        spin_components (numpy array): Contains the excited states spin components ('alp' for alpha spin and 'bet' for beta spin)

    """
    
    critical_params = [ "curr_step" ]
    # Default parameters
    default_params = { "isUKS": 0, "es_software": "cp2k", "logfile_directory": "logfiles" }
    # Check input
    comn.check_input(params, default_params, critical_params)
    
    logfile_directory = params["logfile_directory"]
    es_software = params["es_software"]
    curr_step = params["curr_step"]

    if es_software == "cp2k":
        
        logfile_name = F"{logfile_directory}/step_{curr_step}.log"
        params.update({"logfile_name":logfile_name}) 
        # Extract the excitation energies and their configurations from the output file.
        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = \
        CP2K_methods.read_cp2k_tddfpt_log_file( params )

    elif es_software == "gaussian":

        logfile_name = F"{logfile_directory}/step_{curr_step}.log"
        params.update({"logfile_name":logfile_name}) 
        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = \
        Gaussian_methods.read_gaussian_tddft_log_file( params )
  
    return excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components




def integrate_cube_set( cubefiles_set_1, cubefiles_set_2, dv ):
    """
    This function comutes the overlap matrix between two set of cube files.
    
    Args:
    
        cubefiles_set_1 (list): The list of cube files of the curr_step.

        cubefiles_set_2 (list): The list of cube files of the previous step.

        dv (float): The integration element obtained from the cube files.

    Returns:
    
        overlap_matrix (numpy 2D array): The overlap between the two set of cube files.

    """
    # Initialize the overlap matrix
    # Note: The overlap matrix is a square matrix so the cubefiles_set_1 and cubefiles_set_2 must have the same length
    overlap_matrix = np.zeros( ( len( cubefiles_set_2 ) , len( cubefiles_set_1 ) ) )

    # The overlaps between the cube files of the first set and second set
    for i in range( 0, len( cubefiles_set_2 ) ):
        for j in range( 0, len( cubefiles_set_1 ) ):
            # Use the integrate_cube to compute the integration between the wavefunctions.
            overlap_matrix[i][j] = cube_file_methods.integrate_cube( cubefiles_set_2[i], cubefiles_set_1[j], dv )

    return overlap_matrix





def compute_cube_ks_overlaps( cubefiles_prev, params):
    """
    This function computes overlaps between cube files of two time steps. In order to not read the cube files 
    twice it returns the cube files of the current step, and gets the cube files of the previous step.
    
    Args:
    
        cubefiles_prev (list): The list containing th cube files of the previous step.

        params (dict):

            curr_step (int): The current time step.

            isUKS (int): This parameter is set for spin restricted and unrestricted calculations. When it is
                         set to 1 it means that unrestricted calculations were set in the input file otherwise 
                         it is restricted.
            
            nprocs (int): The number of processors used to read the cube files and perform the integration.
            
    Returns:
    
        cubefiles_curr (list): The list of the read current step cube files.
        
        S_ks_prev (2D numpy array): The overlap matrix of the wavefunctions for the previous time step.
        
        S_ks_curr (2D numpy array): The overlap matrix of the wavefunctions for the current time step.
        
        St_ks (2D numpy array): The overlap matrix between wavefunctions of the two time step.
        
    """
    
    critical_params = [ "curr_step", "nprocs" ]
    # Default parameters
    default_params = { "isUKS": 0, "es_software": "cp2k"}
    # Check input
    comn.check_input(params, default_params, critical_params)
    # Extract the variables
    curr_step = int(params["curr_step"])
    isUKS  = int(params["isUKS"])
    nprocs = int(params["nprocs"])
    es_software = params["es_software"]
    
    print('---------------------------------------------------------')
    print('Starting the calculations by reading the cube files    \n')
    print('                    Step time %i\n'%(curr_step))
    print('-------------------------------------------------------\n')
    print('Reading cube files with %i number of processors'%nprocs)
    # First read all the .cube files (This is the most time consuming part)
    # Creating a pool and assigning the number of processors for mp.Pool
    pool = mp.Pool(processes=nprocs)
    
    # for curr_step
    print("Reading the cubes for step ", curr_step)
    # generate the cube file names produced by CP2K
    cubefile_names_curr = CP2K_methods.cube_file_names_cp2k( params )

    #for cubefile in cubefile_names_curr:
    #    os.system( "/gpfs/scratch/brendan/cp2k/tools/cubecruncher/cubecruncher.x -center geo -i %s -o %s-1.cube " % ( cubefile, cubefile.replace( ".cube", "" ) ) )
    #    os.system( "rm %s" % cubefile)
    #    os.system( "mv %s-1.cube %s" % ( cubefile.replace(".cube",""), cubefile ) )

    # Reading the cube files
    cubefiles_curr = []
    # Apply pool.map to the cube_file_methods to the set of variables of the cubefile_names_curr
    cubefiles_curr = pool.map( cube_file_methods.read_cube, cubefile_names_curr )
    # Close the pool
    pool.close()
    
    # Calculate the dv element for integration
    dv = cube_file_methods.grid_volume( cubefile_names_curr[0] )

    #### Compute the overlaps
    # If spin unrestricted calculations were set to run
    if isUKS == 1:

        # The cube files for alpha and beta spin for the previous time step

        # The alpha cubes are the even indices of the read cube files
        alp_cubes_prev = cubefiles_prev[0::2]
        # The beta cubes are the odd indices of the read cube files
        bet_cubes_prev = cubefiles_prev[1::2]
    
        # The cube files for alpha and beta spin for the current time step

        # The alpha cubes are the even indices of the read cube files
        alp_cubes_curr = cubefiles_curr[0::2]
        # The beta cubes are the odd indices of the read cube files
        bet_cubes_curr = cubefiles_curr[1::2]

        # Initializing the overlap matrices for alpha and beta spin with zero matrices
        zero_mat_alp = np.zeros( ( len( alp_cubes_prev ), len( alp_cubes_curr ) ) )
        zero_mat_bet = np.zeros( ( len( bet_cubes_prev ), len( bet_cubes_curr ) ) )


        ### Using concurrency to compute the integration and overlap matrices
        with concurrent.futures.ThreadPoolExecutor(max_workers=nprocs) as executor:
            # <psi_alp(t-1)|psi_alp(t)>
            int_1 = executor.submit(integrate_cube_set, alp_cubes_prev, alp_cubes_curr, dv )
            # <psi_bet(t-1)|psi_bet(t)>
            int_2 = executor.submit(integrate_cube_set, bet_cubes_prev, bet_cubes_curr, dv )
            # <psi_alp(t-1)|psi_alp(t-1)>
            int_3 = executor.submit(integrate_cube_set, alp_cubes_prev, alp_cubes_prev, dv )
            # <psi_bet(t-1)|psi_bet(t-1)>
            int_4 = executor.submit(integrate_cube_set, bet_cubes_prev, bet_cubes_prev, dv )
            # <psi_alp(t)|psi_alp(t)>
            int_5 = executor.submit(integrate_cube_set, alp_cubes_curr, alp_cubes_curr, dv )
            # <psi_bet(t)|psi_bet(t)>
            int_6 = executor.submit(integrate_cube_set, bet_cubes_curr, bet_cubes_curr, dv )

        # Extracting the results for each of the submitted jobs
        St_alp_alp = int_1.result()
        St_bet_bet = int_2.result()

        # The overlap between cube files at times t-1 and t
        S_alp_alp_prev = int_3.result()
        S_bet_bet_prev = int_4.result()
        S_alp_alp_curr = int_5.result()
        S_bet_bet_curr = int_6.result()

    else:

        ### Using the concurrency
        with concurrent.futures.ThreadPoolExecutor(max_workers=nprocs) as executor:
            # <psi(t-1)|psi(t-1)>
            int_1 = executor.submit(integrate_cube_set, cubefiles_prev, cubefiles_prev, dv )
            # <psi(t)|psi(t)>
            int_2 = executor.submit(integrate_cube_set, cubefiles_curr, cubefiles_curr, dv )
            # <psi(t-1)|psi(t)>
            int_3 = executor.submit(integrate_cube_set, cubefiles_prev, cubefiles_curr, dv )

        # Extracting the results
        S_prev = int_1.result()
        S_curr = int_2.result()
        St     = int_3.result()

        # These are used to form the block matrices for two-spinor format 
        zero_mat_alp = np.zeros( ( len( S_prev ), len( S_curr ) ) )
        zero_mat_bet = np.zeros( ( len( S_prev ), len( S_curr ) ) )

        # The alpha and beta spin have the same overlap matrix in spin restricted case
        S_alp_alp_prev, S_bet_bet_prev = S_prev, S_prev
        S_alp_alp_curr, S_bet_bet_curr = S_curr, S_curr
        St_alp_alp, St_bet_bet = St, St
    
    # Storing the overlap matices in two-spinor format by forming the block matrices
    S_ks_prev = data_conv.form_block_matrix( S_alp_alp_prev, zero_mat_alp, zero_mat_bet, S_bet_bet_prev )
    S_ks_curr = data_conv.form_block_matrix( S_alp_alp_curr, zero_mat_alp, zero_mat_bet, S_bet_bet_curr )
    St_ks     = data_conv.form_block_matrix( St_alp_alp,     zero_mat_alp, zero_mat_bet, St_bet_bet )

    return cubefiles_curr, S_ks_prev, S_ks_curr, St_ks 






def reindex_cp2k_sd_states( ks_orbital_homo_index, ks_orbital_indicies, sd_basis_states ):
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

    ks_orbs_new_index = []
    for i in range(n_alp_ks_orbs):
        ks_orbs_new_index.append(i+1)

    print("\nWe are about to form the excited state SDs, printing the sd_basis_states")
    print(sd_basis_states)
   
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
        sd_excitation = []
        # sd_excitation = [ excitations[j][1] if x == excitations[j][0] else x for x in sd_basis[0] ]
        for sd_state in sd_basis[0]:
            if sd_state==excitations[j][0]:
                sd_excitation.append(excitations[j][1])
            else:
                sd_excitation.append(sd_state)

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



def form_Hvib_real( params ):
    """
    This function forms the real part of the vibronic Hamiltonian by inserting the 
    energies on the diagonal of a zero matrix. The Hvib is in two-spinor format 
    containing the alpha and beta states energies.
    
    Args:
        
        params (distionary):

            logfile_directory (str): The directory of the output files by CP2K.

            time (int): The time step of the molecular dynamics. For single point calculations
                        it should be set to 0.
                    
            min_band (int): The minimum state number.

            max_band (int): The maximum state number.
 
            isUKS (int): This flag is used whenever the UKS calculations were called in CP2K input. if it is set to
                         1 then UKS calculations were called and if not it will consider only alpha energies.
        
    Returns:
    
        Hvib_ks_re (2D numpy array): The diagonal matrix containing the energies of alpha and beta spins if
                                     UKS is set to True and the energies of alpha spin in a block form matrix.
        
    """
    critical_params = [ "min_band", "max_band", "curr_step" ]
    # Default parameters
    default_params = { "isUKS": 0, "es_software": "cp2k", "logfile_directory": "logfiles"}
    # Check input
    comn.check_input(params, default_params, critical_params)
    # Extracting the data

    isUKS = int( params["isUKS"] )
    # ks_orbital_indicies = params["ks_orbital_indicies"]
    es_software = params["es_software"]
    # minimum state
    min_band = params["min_band"] # ks_orbital_indicies[0]
    # maximum state
    max_band = params["max_band"] # ks_orbital_indicies[-1]
    # log file directory
    logfile_directory = params["logfile_directory"]
    # current time step
    curr_step = params["curr_step"]
    # generate the log file name
    logfile_name = F"{logfile_directory}/step_{curr_step}.log"
    # update the logfile_name parameter in params
    params.update({"logfile_name":logfile_name})

    # If the UKS calculations were set to True in CP2K input
    if isUKS == 1:
        # Read energies with alpha spin
        if es_software == "cp2k":
            params["spin"] = 1
            E_alpha, total_energy = CP2K_methods.read_energies_from_cp2k_log_file( params )
            # Read energies with beta spin
            params["spin"] = 2
            E_beta, total_energy = CP2K_methods.read_energies_from_cp2k_log_file( params )
        if es_software == "gaussian":
            params["spin"] = 1
            E_alpha, total_energy = Gaussian_methods.read_energies_from_gaussian_log_file( params )
            # Read energies with beta spin
            params["spin"] = 2
            E_beta, total_energy = Gaussian_methods.read_energies_from_gaussian_log_file( params )
        # Now forming the diagonal matrix containing the Kohn-Sham energies of the alpha and beta spins
        Hvib_ks_re = np.diag( np.concatenate( ( E_alpha, E_beta ) ) )
    # If there is no UKS calculations set
    else:
        if es_software == "cp2k":
            params["spin"] = 1
            # spin = 1
            E_alpha, total_energy = CP2K_methods.read_energies_from_cp2k_log_file( params )
        elif es_software == "gaussian":
            params["spin"] = 1
            E_alpha, total_energy = Gaussian_methods.read_energies_from_gaussian_log_file( params )
        Hvib_ks_re = np.diag( np.concatenate( ( E_alpha, E_alpha ) ) )
    
    return Hvib_ks_re, total_energy



def run_step2_many_body( params ):
    """
    This function is the main function which runs the following calculations:

    * performin CP2K calculations
    * computing the overlap matrices
    * molecular orbital energies
    * nonadiabatic couplings in Kohn-Sham and TD-DFPT level of theory
    * molecular orbital visualization using VMD
    * phase correction and state reordering

    then it will print out all the outputs.

    * Note: The current code works with CP2K but it can be generalized to other electronic structure
            calculation software packages which are capable of producing cube files and performing 
            TD-DFT calculations such as Gaussian. So this parameter is for future implementations.

    Args:

        params (dictionary):

            min_band (integer): The minimum state number.

            max_band (integer): The maximum state number.

            ks_orbital_homo_index (integer): The index of the Kohn-Sham HOMO orbital (State numbers start from 1).

            nsteps_this_job (integer): The number of steps in this job.

            es_software (string): The electronic structure calculation software.

                    * Note: The current code works with CP2K but it can be generalized to other electronic structure
                            calculation software packages which are capable of producing cube files and performing 
                            TD-DFT calculations such as Gaussian. So this parameter is for future implementations.

            es_software_exe (string): The CP2K executable file or path.

            es_software_input_template (str): The full path to CP2K input template.

            project_name (string): The project name for CP2K input.

            isUKS (integer): This is the flag for unrestricted spin calculations. When it is set to
                             1 it will perform the nrestricted calculations otherwise it will consider the restricted case.

            trajectory_xyz_filename (string): The full path to trajectory xyz file.

            njob (integer): The current job number.

            nprocs (integer): The number of processors to be used.

            res_dir (string): The full path to res directory where the output data are stored.

            dt (float): The time step in atomic unit used in molecular dynamics calculations.

            logfile_directory (string): The path to where the log files are stored.

            number_of_states (integer): The number of excited states to be considered in the TD-DFPT calculations.

            tolerance (float): This float number is used to consider only the states which their 
                               square CI coefficients is larger than tolerance factor.

            do_phase_corrections (integer): This flag is for phase correction of the overlap matrices. If it is set to
                                            1 it will perform phase correction otherwise it will not.

            

    Returns:

        None

    """
    # Current working directory
    pwd = os.getcwd()
    # Critical variables
    critical_params = ["min_band", "max_band", "ks_orbital_homo_index", "nsteps_this_job", 'trajectory_xyz_filename']
    # Default parameters
    default_params = { "isUKS": 0, "es_software": "cp2k", "es_software_exe": "cp2k.psmp", "es_software_input_template": "cp2k_input_template.inp", "project_name": "Libra_CP2K", "njob": 1, "nprocs": 2, "dt":  41.3393964448119, "logfile_directory": "logfiles", "number_of_states": 5, "tolerance": 0.0, "do_phase_corrections": 1, "istep": 0, "do_cube_visualization": 0, "states_to_be_plotted": [] }
    # Check input
    comn.check_input(params, default_params, critical_params)
   

    # Extracting the min_band and max_band, we have to use int(...) since 
    # the numbers which are read from the bash are strings.
    min_band = int( params["min_band"] )
    max_band = int( params["max_band"] )
    # Update params with min_band and max_band
    params.update({"min_band":min_band})
    params.update({"max_band":max_band})

    # Define ks_orbital_indicies based on the min_band and max_band    
    ks_orbital_indicies = list( range( min_band, max_band+1 ) )
    # Update params with ks_orbital_indicies
    params.update({"ks_orbital_indicies":ks_orbital_indicies})

    # The HOMO orbital index 
    ks_orbital_homo_index = int(params["ks_orbital_homo_index"])

    # Number of steps for this job
    nsteps_this_job = int(params["nsteps_this_job"])
    es_software = params["es_software"].lower()
    params.update( {"es_software": es_software} )

    # The electronic structure software, for example cp2k
    if es_software == "cp2k":
        # The path to executable CP2K
        cp2k_exe = params["es_software_exe"]
        cp2k_input_template = params["es_software_input_template"]

    elif es_software == "gaussian":
        gaussian_exe = params["es_software_exe"]
        gaussian_input_template = params["es_software_input_template"]

    # The project name, this is important since the cube files are produced based on the project_name
    project_name = params["project_name"]
    # The path to the pre-computed trajectory
    trajectory_xyz_filename = params["trajectory_xyz_filename"]
    # The job number
    njob    = int(params["njob"])
    # The number of processors
    nprocs  = int(params["nprocs"])
    # Update params with nprocs
    params.update({"nprocs":nprocs})

    # The path to res directory where the overlap matrices, energies, and NACs are stored
    res_dir = params["res_dir"]

    # The time step in atomic unit, used to compute the NACs
    dt = float(params["dt"])

    # The path to log files produced by CP2K
    logfile_directory = params["logfile_directory"]

    # Number of excited states
    number_of_states = int(params["number_of_states"])
    # The tolerance factor for CI coefficients
    tolerance = float(params["tolerance"])    

    # Phase corection flag
    do_phase_corrections = int(params["do_phase_corrections"])

    # The current time step which here is set to istep
    curr_step = int(params["istep"])
    # Update the curr_step in params
    params.update({ "curr_step": curr_step })

    # Cube file visualization flag
    do_cube_visualization = int(params["do_cube_visualization"])

    # Flag for the completion level
    completion_level = int(params["completion_level"]) 


    #####################################################################################################
    ######################## Initializing lists for storing the data for each job #######################
    #####################################################################################################
    # For empty lists for the KS, SD, and CI bases. These are to be used for collecting data at each step

    # Overlap of the wavefunctions at the same time step in this job
    S_ks_job  = []
    # Overlap of the wavefunctions for two consecutive time steps in this job
    St_ks_job = []
    # Energy levels of the molecular orbitals for each time step in this job
    E_ks_job  = []
    # The total energies for each time step in this job
    total_energies_job = [] 

    # The overlap matrices for Slater determinant basis for each time step in this job
    S_sd_job  = []
    St_sd_job = []
    # The unique basis of Slater determinants basis
    sd_basis_states_unique = []

    # The overlap matrices for excited states for each time step in this job
    S_ci_job  = []
    St_ci_job = []
    # The configuration interaction coefficients, energies, and basis states with their spin components for each time step in this job
    ci_coefficients_job   = []
    ci_basis_states_job   = []
    ci_energies_job       = []
    spin_components_job   = []

    ####### Setting up the calculations for the initial step

    timer1 = time.time()
    # Set up an initial parameter to compute the norm of the vector

    if es_software == "cp2k":

        # Set up the cp2k input template with the atomic positions for the first timestep of this job batch
        CP2K_methods.CP2K_input_static( cp2k_input_template, project_name, trajectory_xyz_filename, curr_step )
        # Running the CP2K calculations
        os.system("mpirun -np %d %s -i %s-%i.inp -o logfiles/step_%d.log "%( nprocs, cp2k_exe, project_name, curr_step, curr_step ) )
        # After finishing the CP2K calculations move all the pdos and cube files to specified folders
        os.system("mv *.pdos pdosfiles")

    elif es_software == "gaussian":

        Gaussian_methods.gaussian_input( project_name, curr_step, gaussian_input_template, trajectory_xyz_filename )
        os.system("%s < %s-%i.gjf > logfiles/step_%d.log "%( gaussian_exe, project_name, curr_step, curr_step ) )
        Gaussian_methods.cube_generator_gaussian(project_name,curr_step,ks_orbital_indicies[0],ks_orbital_indicies[-1],nprocs,'../../sample_cube_file.cube',int(params["isUKS"]))

    os.system("mv *.cube cubefiles")
    # Print the elapsed time for CP2K calculations for this step.
    print(params["es_software"]," calculation time for step ", params["curr_step"]," was ", time.time() - timer1)

    # We have compted the first SCF calculation for this job, now to read the output data and cubes
    print("Reading the initial step using pool")
    # Creating a pool with nprocs
    pool = mp.Pool( processes = nprocs )

    # The cube file names produced by CP2K, Here we set it as prev since we
    # don't want to read the cubes twice.
    cubefile_names_prev = CP2K_methods.cube_file_names_cp2k( params )
    # We may have to "crunch" the cubes - this may be especially needed for periodic systems.
    #for cubefile in cubefile_names_prev:
    #    os.system( "/gpfs/scratch/brendan/cp2k/tools/cubecruncher/cubecruncher.x -center geo -i %s -o %s-1.cube " % ( cubefile, cubefile.replace( ".cube", "" ) ) )
    #    os.system( "rm %s" % cubefile)
    #    os.system( "mv %s-1.cube %s" % ( cubefile.replace(".cube",""), cubefile ) )
    cubefiles_prev = pool.map( cube_file_methods.read_cube, cubefile_names_prev )

    print("Initial step for this job, the cube file names are:\n", cubefile_names_prev)
    # Call the pool to read the cube files
    # Close the pool 
    pool.close()

    # If we are to visualize the cubefiles, then we should compute the grid-mesh before hand to save time later
    if do_cube_visualization == 1:
        # This phase_factor_visual contains the phase factor correction for each state
        # and is used to plot the cubes phase corrected
        phase_factor_visual = np.ones( ( max_band - min_band + 1 ) * 2 )
        # Update the params with phase_factor_visual
        params.update({"phase_factor_visual":phase_factor_visual})
        # The states which were set to be plotted
        states_to_be_plotted = []
        for state in list( params["states_to_be_plotted"].split(",") ):
            states_to_be_plotted.append( int( state ) )
        # Update the params for states_to_be_plotted
        params.update({"states_to_be_plotted":states_to_be_plotted})
        # Plot the cubes using VMD
        cube_file_methods.plot_cubes( params )

    # After reading the cubefiles for curr_step, we should delete them to be memory efficient
    os.system("rm cubefiles/*")

    # Froming the hvib_ks_re with energies for the first step in the job.
    # This will extract the Kohn-Sham energies and total energy from the CP2K log files and forms the Hvib_real
    hvib_ks_re, total_energy = form_Hvib_real( params )
    # Extracting the energies from Hvib_real matrices
    E_ks_re = data_conv.nparray2CMATRIX( hvib_ks_re )
    # Appending the Kohn-Sham and total energies of this step in the E_ks_job and total_energies_job
    E_ks_job.append( E_ks_re )
    total_energies_job.append( total_energy )

    if completion_level > 0:
        # Get the excitation analysis output
        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = get_excitation_analysis_output( params )
        # Normalize CI coefficients
        ci_coefficients_raw_norm = normalize_ci_coefficients(ci_coefficients_raw_unnorm)
        print("\nPrinting excitation_analysis_output")
        print("excitation_energies = ", excitation_energies)
        print("ci_basis_raw = ", ci_basis_raw)
        print("ci_coefficients_raw_unnorm = ", ci_coefficients_raw_unnorm)
        print("spin_components = ", spin_components)
        #sys.exit(0)

        # Extract the uniquie SD basis states from the ci basis states
        for ci_basis_state_index in range( len( ci_basis_raw ) ):
            for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):
                sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ]
                if sd_basis_state_and_spin not in sd_basis_states_unique:
                    sd_basis_states_unique.append( sd_basis_state_and_spin )
        # Printing the unique Slater determinants basis - this will help the user to 
        # see which states are considered and check them
        print( "Slater determinant basis states = ", sd_basis_states_unique )

        # Appending the extracted data from excitation analysis in the _job variables
        ci_basis_states_job.append( ci_basis_raw )
        ci_coefficients_job.append( ci_coefficients_raw_norm )
        ci_energies_job.append( excitation_energies )
        spin_components_job.append( spin_components )

    curr_step += 1


    #########################################
    # All other steps after initial step for this job
    for step in range( nsteps_this_job-1 ):

        # Update the params with curr_step
        params.update({ "curr_step":curr_step })

        # Setting up the timer
        timer1 = time.time()

        if es_software == "cp2k":

            # Creating the CP2K input
            CP2K_methods.CP2K_input_static( cp2k_input_template, project_name, trajectory_xyz_filename, curr_step )
            # Run CP2K
            os.system("mpirun -np %d %s -i %s-%i.inp -o logfiles/step_%d.log "%( nprocs, cp2k_exe, project_name, curr_step, curr_step ) )
            # Move all the pdos and cube files to one folder
            os.system("mv *.pdos pdosfiles")

        elif es_software == "gaussian":
            Gaussian_methods.gaussian_input( project_name, curr_step, gaussian_input_template, trajectory_xyz_filename )
            os.system("%s < %s-%i.gjf > logfiles/step_%d.log "%( gaussian_exe, project_name, curr_step, curr_step ) )
            Gaussian_methods.cube_generator_gaussian(project_name,curr_step,ks_orbital_indicies[0],ks_orbital_indicies[-1],nprocs,'../../sample_cube_file.cube',int(params["isUKS"]))

        os.system("mv *.cube cubefiles")

        # Print the CP2K timing
        print("CP2K elapsed time for step ", params["curr_step"]," was ", time.time() - timer1)

        # Forming the hvib_ks_re the same as above
        hvib_ks_re, total_energy = form_Hvib_real( params )
        # Extracting the energies and appending them in _job variables
        E_ks_re = data_conv.nparray2CMATRIX( hvib_ks_re )
        E_ks_job.append( E_ks_re )
        total_energies_job.append( total_energy )

        #============================== Computation of the overlap matrices=============================
        # Now, read in the cube files for steps step-1 and step
        # and form the S_ks and St_ks overlap matricies
        cubefiles_curr, S_ks_prev, S_ks_curr, St_ks = compute_cube_ks_overlaps( cubefiles_prev, params )
        # replace the cubefiles_prev with the cubefiles_curr 
        # with this we do not need to read the cubes twice
        cubefiles_prev = cubefiles_curr
        #===============================================================================================

        # Plot the cube files using VMD
        if do_cube_visualization == 1:
            for row_index in range(len(St_ks)):
                # Multiplying the phase factor for each states which needs to be plotted
                if St_ks[row_index][row_index]<0:
                    phase_factor_visual[row_index] = phase_factor_visual[row_index] * (-1)
            # Update the params with phase_factor_visual
            params.update({"phase_factor_visual":phase_factor_visual})
            # Plot the cubes
            cube_file_methods.plot_cubes( params )

        # Now remove all the cube files to be memory efficient
        os.system("rm cubefiles/*")

        # Appending the normalized overlap matrices in the _job variables
        S_ks_job.append(S_ks_prev)
        St_ks_job.append(St_ks)
        if step == nsteps_this_job-2:
            S_ks_job.append(S_ks_curr)

        if completion_level > 0:
            print("Reading the excitation analysis output for step ", params["curr_step"])
            # Get excitation analysis output results
            excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = get_excitation_analysis_output( params )
            print("\nPrinting excitation_analysis_output")
            print("excitation_energies = ", excitation_energies)
            print("ci_basis_raw = ", ci_basis_raw)
            print("ci_coefficients_raw_unnorm = ", ci_coefficients_raw_unnorm)
            print("spin_components = ", spin_components)
            # Normalize the coefficients
            ci_coefficients_raw_norm = normalize_ci_coefficients(ci_coefficients_raw_unnorm)
            print("ci_coefficients_raw_norm = ", ci_coefficients_raw_norm)
            # Extract the uniquie SD basis states from the ci basis states

            for ci_basis_state_index in range( len( ci_basis_raw ) ):
                for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):
                    sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ] 
                    if sd_basis_state_and_spin not in sd_basis_states_unique:
                        sd_basis_states_unique.append( sd_basis_state_and_spin )
            print( "The unique SD basis states = ", sd_basis_states_unique )

            # Now append the extracted excitation analysis output in _job variables
            ci_basis_states_job.append( ci_basis_raw )
            ci_coefficients_job.append(   ci_coefficients_raw_norm )
            ci_energies_job.append( excitation_energies )
            spin_components_job.append( spin_components )

        curr_step += 1

    # We need to output the last step overlap and energies as well.
    step = nsteps_this_job-1
    for step in range( nsteps_this_job ):
        S_ks_job[step]  = data_conv.nparray2CMATRIX(  S_ks_job[step] )
    for step in range( nsteps_this_job-1 ):
        St_ks_job[step] = data_conv.nparray2CMATRIX( St_ks_job[step] )
    step3.apply_normalization( S_ks_job, St_ks_job )
    for step in range( nsteps_this_job ):
        S_ks_job[step].real().show_matrix("%s/S_ks_%d_re" % (res_dir, int(params["istep"])+step))
        E_ks_job[step].real().show_matrix("%s/E_ks_%d_re" % (res_dir, int(params["istep"])+step))
    for step in range( nsteps_this_job-1 ):
        St_ks_job[step].real().show_matrix("%s/St_ks_%d_re" % (res_dir, int(params["istep"])+step))
    #sys.exit(0)

    if completion_level == 0:
        print("\nComplete! Exiting for completion levels 0 or 1")
        sys.exit(0)

    #################################################################################################################
    #################################################################################################################
    print("Finished with all of the step. Now computing the overlap matrices and NACs in TD-DFPT level of theory...")
    # Now, time to compute S_sd and St_sd
    # Start by reindexing the unique Slater determinant basis. The current SD bases are not able to be read by Libra
    sd_states_reindexed = reindex_cp2k_sd_states( ks_orbital_homo_index, ks_orbital_indicies, sd_basis_states_unique )
    print("The reindexed Slater determinants = \n", sd_states_reindexed)

    # For each step, we must sort the set of unique Slater determinant states based on their energy
    E_sd_job = []
    sd_states_reindexed_sorted = []
    sd_states_unique_sorted = []
    SD_energy_corr = [0.0]*len(sd_states_reindexed)
    for step in range( nsteps_this_job ):
        E_this_sd  = mapping.energy_mat_arb( sd_states_reindexed, E_ks_job[step], SD_energy_corr )
        nstates_sd = len(sd_states_reindexed)
        e = np.zeros( nstates_sd )
        for state in range(nstates_sd):
            e[state] =  E_this_sd.get(state,state).real
        reindex = np.argsort(e)
        E_sd_job.append(  CMATRIX(nstates_sd,nstates_sd) )
        sd_states_reindexed_sorted.append( [] )
        for i in range(len(reindex)):
            E_sd_job[step].set(  i,i, E_this_sd.get(  int(reindex[i]), int(reindex[i])) )
            sd_states_reindexed_sorted[step].append( sd_states_reindexed[ int(reindex[i]) ] )
        # To omit the ground state, we will manually make this later. the below variable goes to make the T matrix
        sd_states_unique_sorted.append( [] )
        for i in range(1,len(reindex)):
            sd_states_unique_sorted[step].append( sd_basis_states_unique[ int(reindex[i])-1 ] )

    # For each step make S_sd and St_sd
    print("Making the S_sd and St_sd....")
    for step in range( nsteps_this_job ):
        s_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step], S_ks_job[step],  use_minimal=False )
        S_sd_job.append( s_sd )
    for step in range( nsteps_this_job-1 ):
        st_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step+1], St_ks_job[step],  use_minimal=False )
        St_sd_job.append( st_sd )

    # Output Slater determinant data to the res directory
    print("Outputting the SD overlaps and energies.")
    for step in range( nsteps_this_job ):
        S_sd_job[step].real().show_matrix("%s/S_sd_%d_re" % (res_dir, int(params["istep"])+step))
        E_sd_job[step].real().show_matrix("%s/E_sd_%d_re" % (res_dir, int(params["istep"])+step))
    for step in range( nsteps_this_job-1 ):
        St_sd_job[step].real().show_matrix("%s/St_sd_%d_re" % (res_dir, int(params["istep"])+step))

    # Now, we have computed the overlaps in the KS and SD bases. Now, we need to compute the SD2CI matrix for each step
    # Now, form the SD2CI matrix for each step
    ci_coefficients = []
    # Add one to the number of CI states because the ground state is not included yet
    nSDs = len( sd_basis_states_unique ) + 1
    nCIs = number_of_states + 1
    SD2CI = []

    print("Computing the NACs for single-particle and many-body wavefunctions....")
    for step in range( nsteps_this_job ):

        # Make the list of ci_coefficients for each step in the way Libra accepts
        ci_coefficients.append( [] )
        # Start with the ground state. This is not explicitly given by electronic strcture calculations
        ci_coefficients[step].insert( 0, [0.0] * nSDs )
        ci_coefficients[step][0][0] = 1.0

        # For each ci state for this step
        for i in range( len( ci_coefficients_job[step] ) ):
            count = 0
            ci_coefficients[step].append( [0.0] * nSDs )
            # Exclude ground state here in the index, that info is not explicitly contained 
            # in the ci_coefficients_dynamics list from electronic strcture calculations
            tmp_ci_basis_state_and_spin = []
            for k in range(len(ci_coefficients_job[step][i])):
                tmp_ci_basis_state_and_spin.append( [ci_basis_states_job[step][i][k] , spin_components_job[step][i][k]] )

            for j in range( nSDs-1 ):
                if sd_states_unique_sorted[step][j] in tmp_ci_basis_state_and_spin:   
                    # ok, it has found a match, now what is the index?
                    item_index = tmp_ci_basis_state_and_spin.index(sd_states_unique_sorted[step][j])
                    ci_coefficients[step][i+1][j+1] = float(ci_coefficients_job[step][i][item_index])    
 
        SD2CI.append( CMATRIX( nSDs, nCIs ) )
        for i in range( nSDs ):
            for j in range( nCIs ):
                SD2CI[step].set( i, j, ci_coefficients[step][j][i] * (1.0+0.0j) )
        SD2CI[step].show_matrix( "T_%s.txt" % str(step) )

    # For each step make S_ci and St_ci
    print("Making the S_ci and St_ci matrices....")
    for step in range( nsteps_this_job ):
        s_ci = SD2CI[step].H()  * S_sd_job[step]  * SD2CI[step]
        S_ci_job.append( s_ci )
    for step in range( nsteps_this_job-1 ):
        st_ci = SD2CI[step].H() * St_sd_job[step] * SD2CI[step+1]
        St_ci_job.append( st_ci )

    # Now, compute the CI energy matrix at each-point and the mid-points
    # For each step
    print("Computing the CI energy matrices....")
    ci_energies_job_cmatrix = []
    for step in range( nsteps_this_job ):
        ci_energies_job_cmatrix.append( CMATRIX( number_of_states + 1, number_of_states + 1 ) )
        for state in range( number_of_states + 1 ):
            if state == 0:
                ci_energies_job_cmatrix[step].set( state, state, total_energies_job[step] )
            else:
                ci_energies_job_cmatrix[step].set( state, state, total_energies_job[step] + ( ci_energies_job[step][state-1]  * units.ev2Ha )  )

    # At the midpoints
    ci_midpoint_energies = []
    for step in range( nsteps_this_job-1 ):
        total_energy_mid_point = 0.5 * ( total_energies_job[step] + total_energies_job[step+1] )
        ci_midpoint_energies.append( CMATRIX( number_of_states + 1, number_of_states + 1 ) )
        for state in range( number_of_states + 1 ):
            if state == 0:
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point )
            else:
                midpoint_energy = 0.5 * ( ci_energies_job[step][state-1] + ci_energies_job[step+1][state-1] )
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point + ( midpoint_energy  * units.ev2Ha )  )

    # Are we to perform state reordering?
    if int(params["perform_state_reordering"]) == 1:
        params2 = {"do_state_reordering":int(params["do_state_reordering"]), "state_reordering_alpha":float(params["state_reordering_alpha"])}
        print("Applying state reordering....")
        step3.apply_state_reordering_general( St_ci_job, ci_midpoint_energies, params2 )

    # Are we to perform phase corrections?
    if do_phase_corrections == 1:
        print("\nApplying phase corrections")
        step3.apply_phase_correction_general( St_ci_job )

    # Output CI data to res directory
    print("Outputting the CI data to the res directory..." )
    for step in range( nsteps_this_job ):
        S_ci_job[step].real().show_matrix("%s/S_ci_%d_re"   % (res_dir, int(params["istep"])+step))
        ci_energies_job_cmatrix[step].real().show_matrix("%s/E_ci_%d_re"   % (res_dir, int(params["istep"])+step))
    for step in range( nsteps_this_job-1 ):
        St_ci_job[step].real().show_matrix("%s/St_ci_%d_re" % (res_dir, int(params["istep"])+step))

    # Now, compute the CI NACs and compute the CI Hvib
    print("Computing and outputting the CI NACs...")
    for step in range( nsteps_this_job-1 ): 
        ci_nacs = -(  0.5j / dt ) * CMATRIX ( ( St_ci_job[step] - St_ci_job[step].H() ).real() )    
        ci_hvib = ci_midpoint_energies[step] + ci_nacs
        ci_hvib.real().show_matrix("%s/Hvib_ci_%d_re" % (res_dir, int(params["istep"])+step))
        ci_hvib.imag().show_matrix("%s/Hvib_ci_%d_im" % (res_dir, int(params["istep"])+step))


    print("All steps were done successfully for this job!")

 
