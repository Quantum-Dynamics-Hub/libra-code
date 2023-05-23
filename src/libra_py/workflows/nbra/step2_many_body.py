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
from distutils.spawn import find_executable
import multiprocessing as mp

# Concurrency for large systems
import concurrent.futures

from liblibra_core import *

from libra_py import data_conv
from libra_py import cube_file_methods
import libra_py.packages.cp2k.methods as CP2K_methods
import libra_py.packages.gaussian.methods as Gaussian_methods
import libra_py.packages.dftbplus.methods as DFTB_methods
from libra_py.workflows.nbra import mapping
from libra_py.workflows.nbra import step3
from libra_py import units

import util.libutil as comn


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

    elif es_software == "dftb+":

        logfile_name = F"{logfile_directory}/TRA_{curr_step}.DAT"
        params.update({"logfile_name":logfile_name})
        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = \
        DFTB_methods.read_dftbplus_TRA_file( params )
  
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




def reindex_cp2k_sd_states( ks_orbital_homo_index, ks_orbital_indicies, sd_basis_states, sd_format=2 ):
    """
    ks_orbital_homo_index: Index of the homo ks orbital, from 1
    ks_orbital_indicies: Range of the considered ks orbtials. Ex) [8,9,10,11], where 9 is homo orbtial index (from 1)
    sd_basis_states( list of lists of lists ): A list of Slater determinants, where each slater determinant is a excitation in the Kohn-Sham
                                               basis. This function assumes that all Kohn-Sham excitations are for alpha electrons. To
                                               differentiate between alpha and beta excitations, elements of sd_basis_states contain spin
                                               information.
               
                                               Ex) sd_basis_states[0] = [ [9,10], "alp" ] 
                                                   sd_basis_states[1] = [ [9,10], "bet" ]
    Returns:
         Reindexed Slater determinant basis states in terms of the ks_orbital_homo_index. The SD states returned can be either in one of two ways

         1. All alpha orbitals before beta orbtials. This format is still a fixed slot format. It is used when used both alpha and beta spin channels
         Ex.)   ks_orbital_homo_index = 68,    ks_orbital_indicies = [67, 68, 69, 70, 71, 72]
                The ground state will be: [1, 2, -7, -8]
                [ [68,69], "alp" ] --> [1, 3, -7, -8]
                [ [68,69], "bet" ] --> [1, 2, -7, -9]

         2. Fixed slot format. this format used only the alpha spin channel
         Ex.)   ks_orbital_homo_index = 68,    ks_orbital_indicies = [67, 68, 69, 70, 71, 72]
                The ground state will be: [1, -1, 2, -2]
                [ [68,69], "alp" ] --> [1, -1, 3, -2]
                [ [68,69], "bet" ] --> [1, -1, 2, -3]

    """

    #print("\nEntered reindex_cp2k_sd_states function")
    #print("ks_orbital_homo_index", ks_orbital_homo_index)
    #print("ks_orbital_indicies", ks_orbital_indicies)
    #print("sd_basis_states", sd_basis_states)

    # We need to update the indexing of the sd_basis - in terms of the rows and cols of St_KS
    # reindex ks orbs according to the matrix size
    n_alp_ks_orbs = len(ks_orbital_indicies)
    alp_homo_matrix_index = 0
    for i in range( n_alp_ks_orbs ):
        if ks_orbital_indicies[i] == ks_orbital_homo_index:
            alp_homo_matrix_index = i

    #print("alp_homo_matrix_index",alp_homo_matrix_index)

    ks_orbs_new_index = []
    for i in range( n_alp_ks_orbs ):
        ks_orbs_new_index.append( i+1 )

    # Form excited state SDs
    excitations = []
    # For each Slater determinant basis state, which could have spin_component "alp" or "bet"
    for j in range( len( sd_basis_states ) ):

        #print(int( sd_basis_states[j][0][0] ))
        #print(int( sd_basis_states[j][0][1] ))

        if sd_format == 1:
            if sd_basis_states[j][1] == "alp":
                initial_ks_orb = int( sd_basis_states[j][0][0] ) - ks_orbital_homo_index + alp_homo_matrix_index + 1
                final_ks_orb   = int( sd_basis_states[j][0][1] ) - ks_orbital_homo_index + alp_homo_matrix_index + 1
            elif sd_basis_states[j][1] == "bet":
                initial_ks_orb = int( sd_basis_states[j][0][0] ) - ks_orbital_homo_index + alp_homo_matrix_index + n_alp_ks_orbs + 1
                final_ks_orb   = int( sd_basis_states[j][0][1] ) - ks_orbital_homo_index + alp_homo_matrix_index + n_alp_ks_orbs + 1

        elif sd_format == 2:  
            initial_ks_orb = int( sd_basis_states[j][0][0] ) - ks_orbital_homo_index + alp_homo_matrix_index + 1
            final_ks_orb   = int( sd_basis_states[j][0][1] ) - ks_orbital_homo_index + alp_homo_matrix_index + 1

        if sd_basis_states[j][1] == "alp":
            excitations.append( [initial_ks_orb, final_ks_orb]  )
    
        elif sd_basis_states[j][1] == "bet":
            excitations.append( [-initial_ks_orb, -final_ks_orb]  )

    #print( "excitations = ", excitations )
    #sys.exit(0)

    # Form ground-state SD first
    sd_basis = [ [] ]
    if sd_format == 1:
        for i in range( 1, 2*len( ks_orbital_indicies ) ):
            if i < alp_homo_matrix_index + 2:
                sd_basis[0].append( i )
            if i > n_alp_ks_orbs and i < n_alp_ks_orbs + alp_homo_matrix_index + 2:
                sd_basis[0].append( -i )

    elif sd_format == 2:
        for i in range( 1, len( ks_orbital_indicies ) ):
            if i < alp_homo_matrix_index + 2:
                sd_basis[0].append(  i )
                sd_basis[0].append( -i )
    #print ( "ground state = ", sd_basis )

    # Now that we have done the ground state slater 
    for j in range( len( excitations ) ):
        sd_excitation = []
        for sd_state in sd_basis[0]:
            if sd_state==excitations[j][0]:
                sd_excitation.append(excitations[j][1])
            else:
                sd_excitation.append(sd_state)
        sd_basis.append( sd_excitation )
    #print ( sd_basis )
    
    return sd_basis




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
        elif es_software == "dftb+":
            # generate the log file name
            logfile_name = F"{logfile_directory}/step_{curr_step}_ks.log"
            # update the logfile_name parameter in params
            params.update({"logfile_name":logfile_name})
            params["spin"] = 1
            E_alpha, total_energy = DFTB_methods.get_dftb_ks_energies( params )
        Hvib_ks_re = np.diag( np.concatenate( ( E_alpha, E_alpha ) ) )
    
    return Hvib_ks_re, total_energy



def run_step2_many_body( params ):
    """
    This function is the main function which runs the following calculations:

    * Executes calls to perform CP2K, DFTB+, or Gaussian electronic structure calculations
    * Computing the overlap, time-overlap, and energy eigenvalue matrices for the KS basis
    * Molecular orbital visualization using VMD

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

            logfile_directory (string): The path to where the log files are stored.

    Returns:

        None

    """
    # Current working directory
    pwd = os.getcwd()
    # Critical variables
    critical_params = []
    # Default parameters
    default_params = { "min_band":1, "max_band":1, "ks_orbital_homo_index":0, 
                       "nsteps_this_job":1, 'trajectory_xyz_filename':"md.xyz", "isUKS": 0, 
                       "es_software": "cp2k", "es_software_exe": "cp2k.popt", 
                       "es_software_input_template": "cp2k_input_template.inp", 
                       "project_name": "Libra_CP2K", "njob": 1, "nprocs": 2, "logfile_directory": "logfiles", 
                       "istep": 0, "do_cube_visualization": 0, "states_to_be_plotted": [], 
                       "waveplot_exe":"/util/academic/dftbplus/20.2.1-arpack/bin/waveplot" }
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

    elif es_software == "dftb+":
        dftbp_exe = params["es_software_exe"]
        dftb_input_template = params["es_software_input_template"]
        waveplot_exe = params["waveplot_exe"]

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

    # The path to log files produced by CP2K
    logfile_directory = params["logfile_directory"]

    # The current time step which here is set to istep
    curr_step = int(params["istep"])
    # Update the curr_step in params
    params.update({ "curr_step": curr_step })

    # Cube file visualization flag
    do_cube_visualization = int(params["do_cube_visualization"])

    # Make a directory for this job folder for storing the logfiles, cubefiles, and pdosfiles
    if not os.path.exists("logfiles"):
        os.mkdir("logfiles")
    if not os.path.exists("cubefiles"):
        os.mkdir("cubefiles")
    if not os.path.exists("pdosfiles"):
        os.mkdir("pdosfiles")
    if not os.path.exists("../../all_logfiles"):
        os.mkdir("../../all_logfiles")
    if not os.path.exists("../../all_pdosfiles"):
        os.mkdir("../../all_pdosfiles")

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

    elif es_software == "dftb+":

        print( 'params["curr_step"]', params["curr_step"])
        # Set up the dftb+ input template.
        DFTB_methods.make_dftb_input( dftb_input_template, 0, params["curr_step"]  )
        # Convert the xyz file for currstep to the .gen file format
        os.system("xyz2gen coord-"+str(params["curr_step"])+".xyz")
        # Run the dftb+ program
        if find_executable('srun') is not None:
            os.system("srun %s"%( dftbp_exe ) )
        else:
            os.system("%s"%( dftbp_exe ) )
        os.system("mv band.out step_"+str(curr_step)+"_ks.log")
        os.system("mv step_"+str(curr_step)+"_ks.log logfiles/.")
        os.system("mv EXC.DAT EXC_"+str(curr_step)+".DAT")
        os.system("mv EXC_"+str(curr_step)+".DAT logfiles/.")
        os.system("mv TRA.DAT TRA_"+str(curr_step)+".DAT")
        os.system("mv TRA_"+str(curr_step)+".DAT logfiles/.")
        DFTB_methods.cube_generator_dftbplus( project_name, curr_step, ks_orbital_indicies[0], ks_orbital_indicies[-1], waveplot_exe, int(params["isUKS"]) )

    os.system("mv *.cube cubefiles")
    # Print the elapsed time for CP2K calculations for this step.
    print(params["es_software"]," calculation time for step ", params["curr_step"]," was ", time.time() - timer1)

    # We have compted the first SCF calculation for this job, now to read the output data and cubes
    print("Reading the initial step using pool")
    print("nprocs", nprocs)
    # Creating a pool with nprocs
    pool = mp.Pool( processes = nprocs )

    # The cube file names produced by CP2K, Here we set it as prev since we
    # don't want to read the cubes twice.
    cubefile_names_prev = CP2K_methods.cube_file_names_cp2k( params )
    print("Initial step for this job, the cube file names are:\n", cubefile_names_prev)
    # We may have to "crunch" the cubes - this may be especially needed for periodic systems.
    #for cubefile in cubefile_names_prev:
    #    os.system( "/gpfs/scratch/brendan/cp2k/tools/cubecruncher/cubecruncher.x -center geo -i %s -o %s-1.cube " % ( cubefile, cubefile.replace( ".cube", "" ) ) )
    #    os.system( "rm %s" % cubefile)
    #    os.system( "mv %s-1.cube %s" % ( cubefile.replace(".cube",""), cubefile ) )
    cubefiles_prev = pool.map( cube_file_methods.read_cube, cubefile_names_prev )

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
    print("Removing cubefiles")
    os.system("rm cubefiles/*")

    # Froming the hvib_ks_re with energies for the first step in the job.
    # This will extract the Kohn-Sham energies and total energy from the CP2K log files and forms the Hvib_real
    hvib_ks_re, total_energy = form_Hvib_real( params )
    # Appending the Kohn-Sham and total energies of this step in the E_ks_job and total_energies_job
    E_ks_job.append( hvib_ks_re )
    total_energies_job.append( total_energy )

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

        elif es_software == "dftb+":

            print( 'params["curr_step"]', params["curr_step"])
            # Set up the dftb+ input template.
            DFTB_methods.make_dftb_input( dftb_input_template, 0, params["curr_step"]  )
            # Convert the xyz file for currstep to the .gen file format
            os.system("xyz2gen coord-"+str(params["curr_step"])+".xyz")
            # Run the dftb+ program
            if find_executable('srun') is not None:
                os.system("srun %s"%( dftbp_exe ) )
            else:
                os.system("%s"%( dftbp_exe ) )
            os.system("mv band.out step_"+str(curr_step)+"_ks.log")
            os.system("mv step_"+str(curr_step)+"_ks.log logfiles/.")
            os.system("mv EXC.DAT EXC_"+str(curr_step)+".DAT")
            os.system("mv EXC_"+str(curr_step)+".DAT logfiles/.")
            os.system("mv TRA.DAT TRA_"+str(curr_step)+".DAT")
            os.system("mv TRA_"+str(curr_step)+".DAT logfiles/.")
            DFTB_methods.cube_generator_dftbplus( project_name, curr_step, ks_orbital_indicies[0], ks_orbital_indicies[-1], waveplot_exe, int(params["isUKS"]) )
  
        os.system("mv *.cube cubefiles")
        # Print the Timing
        print("Elapsed time for step ", params["curr_step"]," was ", time.time() - timer1)

        # Forming the hvib_ks_re the same as above
        hvib_ks_re, total_energy = form_Hvib_real( params )
        E_ks_job.append( hvib_ks_re )
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

        #np.savetxt(filename, np.array, fmt='%.16e', delimiter=" ")
        np.savetxt("%s/S_ks_%d_re"  % (res_dir, int(params["istep"])+step), S_ks_job[step], fmt='%.16e', delimiter=" ")
        np.savetxt("%s/E_ks_%d_re"  % (res_dir, int(params["istep"])+step), E_ks_job[step], fmt='%.16e', delimiter=" ")
        np.savetxt("%s/St_ks_%d_re" % (res_dir, int(params["istep"])+step), St_ks_job[step], fmt='%.16e', delimiter=" ")

        curr_step += 1

    # Print out the KS Overlap and Energy matricies for the last step in this job batch
    np.savetxt("%s/S_ks_%d_re" % (res_dir, int(params["istep"])+step+1), S_ks_job[step+1], fmt='%.16e', delimiter=" ")
    np.savetxt("%s/E_ks_%d_re" % (res_dir, int(params["istep"])+step+1), E_ks_job[step+1], fmt='%.16e', delimiter=" ")
    print("All steps were done successfully for this job!")

    os.system("mv logfiles/* ../../all_logfiles/.")
    os.system("mv pdosfiles/* ../../all_pdosfiles/.")
    os.system("rm "+trajectory_xyz_filename)
    os.system("rm *wfn")
