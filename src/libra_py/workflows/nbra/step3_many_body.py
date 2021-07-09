#*********************************************************************************
#* Copyright (C) 2021 Alexey V. Akimov, Brendan A. Smith, Mohammad Shakiba
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  
#
"""
.. module:: step3_many_body
   :platform: Unix, Windows
   :synopsis: This module is designed to convert the results of QE calculations 
       (KS orbital energies and time-overlaps in the KS basis) to the generic Hvib
       matrices, which account for:

           - state reordering;
           - phase corrections;
           - multi-electron wavefunction (Slater determinants) and spin-adaptation
           - scissor operator corrections to energy levels

.. moduleauthor:: Alexey V. Akimov, Brendan A. Smith, Mohammad Shakiba

"""

import os, sys, time, math, cmath
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
from libra_py import units as units
from libra_py import influence_spectrum as infsp
from libra_py import data_visualize, data_conv, data_read, data_stat, data_outs
from . import mapping, step2_many_body, step3, step4


#######################################################################################
# Helper functions
def get_step2_mb_sp_properties( _params ):
    """
    This function extracts information from es output files. It currently works for
    cp2k, dftb+, gaussian

    Args: 
        params ( dictionary ) control parameters
        
            params["orbital_indices"] (list)   - ks orbitals that were used in step2. [default: [0] ] Ex) [15,16,17, ... , 32,33,34]
            params["logfile_directory"] (string): [default: "all_logfiles"]
            params["es_software"] (string): software used for calculations. Options: "cp2k", "dftb+" [default: "cp2k"]
            params["isUKS"] (int): 0 for spin-unpolarized, 1 for spin-polarized [ default: 0]
            params["number_of_states"] (int): number of ci states to extract [ default: 1]
            params["tolerance"] (float): cutoff for SD contribution [default: 0.01]
            params["isnap"] (int): file index to start reading from [default: 0]
            params["fsnap"] (int): file index to stop  reading from [default: 1]
   
    Returns:
        tuple: (S_sd_job, St_sd_job, sd_basis_states_unique, ci_basis_states_job, ci_coefficients_job, ci_energies_job, spin_components_job)

        Here:     

            S_sd_job  (list) - overlaps of SD at each step
            St_sd_job (list) - time-overlaps of SD at each step
            sd_basis_states_unique (list) - 1 of each of the SP transitions (and its spin) that made up the considered CI states
            ci_basis_states_job (list) - similar to sd_basis_states_unique, but all SP transitions encountered
            ci_coefficients_job (list) - all ci coefficients encountered
            ci_energies_job (list) - excitation energies
            spin_components_job (list) - just the spin components (alpha or beta excitaiton?)

    """

    params = dict(_params)

    critical_params = []
    # Default parameters
    default_params = { "orbital_indices":[0], "logfile_directory":"all_logfiles",
                       "es_software":"cp2k", "isUKS":0, "number_of_states":1, 
                       "tolerance":0.01, "isnap":0, "fsnap":1 }
    # Check input
    comn.check_input(params, default_params, critical_params)  

    # Just a quick fix to input parameter names
    params["ks_orbital_indicies"] = params["orbital_indices"]


    start_time  = params["isnap"]
    finish_time = params["fsnap"]

    curr_step = start_time

    # Define ks_orbital_indicies based on the min_band and max_band    

    # Update params with ks_orbital_indicies
    params["curr_step"] = curr_step
    params["logfile_name"] = F"step_{curr_step}.out"

    #S_sd_job, St_sd_job = [], []
    ci_energies_job, ci_basis_states_job, ci_coefficients_job = [], [], []
    spin_components_job = []
    sd_basis_states_unique = []

    excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = step2_many_body.get_excitation_analysis_output( params )
    ci_coefficients_raw_norm = step2_many_body.normalize_ci_coefficients(ci_coefficients_raw_unnorm)

    #print( excitation_energies )
    #print( ci_basis_raw )
    #print( ci_coefficients_raw_unnorm )
    #print( spin_components )
    #print( ci_coefficients_raw_norm )

    # Now append the extracted excitation analysis output to the respective lists
    ci_basis_states_job.append( ci_basis_raw )
    ci_coefficients_job.append( ci_coefficients_raw_norm )
    ci_energies_job.append( excitation_energies )
    spin_components_job.append( spin_components )

    # Extract the uniquie SD basis states from the ci basis states
    for ci_basis_state_index in range( len( ci_basis_raw ) ):
        for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):

             sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , \
                                         spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ]
             if sd_basis_state_and_spin not in sd_basis_states_unique:
                 sd_basis_states_unique.append( sd_basis_state_and_spin )

    #print( "Slater determinant basis states = ", sd_basis_states_unique )
    curr_step += 1

    # All other steps after initial step for this job
    for step in range( finish_time-start_time-1 ):

        params.update({ "curr_step":curr_step })
        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = step2_many_body.get_excitation_analysis_output( params )
        ci_coefficients_raw_norm = step2_many_body.normalize_ci_coefficients(ci_coefficients_raw_unnorm)

        # Extract the uniquie SD basis states from the ci basis states
        for ci_basis_state_index in range( len( ci_basis_raw ) ):
            for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):

                sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , \
                                            spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ]
                if sd_basis_state_and_spin not in sd_basis_states_unique:
                    sd_basis_states_unique.append( sd_basis_state_and_spin )
        #print( "Slater determinant basis states = ", sd_basis_states_unique )

        # Now append the extracted excitation analysis output in _job variables
        ci_basis_states_job.append( ci_basis_raw )
        ci_coefficients_job.append(   ci_coefficients_raw_norm )
        ci_energies_job.append( excitation_energies )
        spin_components_job.append( spin_components )
        curr_step += 1

    return sd_basis_states_unique, ci_basis_states_job, ci_coefficients_job, ci_energies_job, spin_components_job
    #return S_sd_job, St_sd_job, sd_basis_states_unique, ci_basis_states_job, ci_coefficients_job, ci_energies_job, spin_components_job






def compute_ci_energies_midpoint( ci_energies, _params ):
    """
    This function compute the excitation energies energies at the midpoint from a list of excitation energies at each step. 
    At each step, there are many electronic states. This function takes a list as an input, and is meant to be used 
    in the NBRA workflow calculatiosn where lists may be more convenient than matricies. 

    This funciton is made to be used within the NBRA Libra workflow, where things such as ci_energies have been extracted from TD-DFT calculations. 
    As of 11/30/2020, compatable ES programs include CP2K, DFTB+ and Gaussian.

    Energies are assumed to be energies from TDDFT calculatons. This function gives zero as the ground state total energy

    Args:
        ci_energies (list of lists): energies of the MB states
        num_excited_states (int): number of excited states
        istep (int): step at which to start counting
        fstep (int): stap at which to stop counting

    Returns:
        ci_midpoint_energies (list of CMATRIX): energies in Ha. Ground state energy is set to zero
    """


    params = dict(_params)

    critical_params = []
    # Default parameters
    default_params = { "number_of_states":1, "isnap":0, "fsnap":1 }
    # Check input
    comn.check_input(params, default_params, critical_params)  


    nstates = params["number_of_states"]
    istep   = params["isnap"]
    fstep   = params["fsnap"]


    # Now, compute the CI energy matrix at each-point and the mid-points
    # For each step
    #print("Computing the CI energy matrices....")
    ci_energies_cmatrix = []
    for step in range( fstep - istep ):
        ci_energies_cmatrix.append( CMATRIX( nstates + 1, nstates + 1 ) )

        for state in range( nstates + 1 ):
            if state == 0:
                ci_energies_cmatrix[step].set( state, state, 0.0 )
            else:
                ci_energies_cmatrix[step].set( state, state, ( ci_energies[step][state-1]  * units.ev2Ha )  )


    # At the midpoints
    ci_midpoint_energies = []
    for step in range( fstep - istep - 1 ):
        total_energy_mid_point = 0.0 #0.5 * ( total_energies[step] + total_energies[step+1] )
        ci_midpoint_energies.append( CMATRIX( nstates + 1, nstates + 1 ) )
        for state in range( nstates + 1 ):
            if state == 0:
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point )
            else:
                midpoint_energy = 0.5 * ( ci_energies[step][state-1] + ci_energies[step+1][state-1] )
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point + ( midpoint_energy  * units.ev2Ha )  )

    
    return ci_midpoint_energies




def make_T_matrices( ci_coefficients, ci_basis_states, spin_components, sd_states_unique_sorted, _params):
    """
    This function makes the "T"ransformation matricies that convert between the SD basis to the CI-like (or many-body (MB)) basis.

    This funciton is made to be used within the NBRA Libra workflow, where things such as ci_coefficients, ci_basis_states, spin_components, 
    and sd_states_unique_sorted have been extracted from TD-DFT calculations. As of 11/30/2020, compatable ES programs
    include CP2K, DFTB+ and Gaussian.

    Args:
        ci_coefficients (list of lists of lists): coefficients for the many-body states for each step
        ci_basis_states (list of lists): All SD basis states that comprise the many-body excitations for each step
        spin_components (list of lists): the spin components of the excitation (alpha or beta excitaiton?) for all states and all steps  
        sd_basis_states_unique (list): 1 of each of the SP transitions (and its spin) that made up the considered CI states
        _params (dict): control parameters
            * **_params["number_of_states"]** (int): number of excited MB states [default: 1]
            * **_params["isnap"] (int): step at which to start counting  [ default: 0 ]
            * **_params["fsnap"] (int): step at which to stop counting [ default: 1]
            * **_params["outdir"] (string): output directory for the T matricies [ default: "res"]
            * **_params["verbosity"] (int): level of extra m messages you wanna see [ default: 0 - no debug info ] 

    Returns:
        SD2CI (list of CMATRIX): CMATRIX at each timestep where the rows are SDs and the cols are MB states. The columns contain the coefficients of the MB expansion for each MB state

    """

    params = dict(_params)

    critical_params = []
    # Default parameters
    default_params = { "number_of_states":1, "isnap":0, "fsnap":1, "verbosity":0, "outdir":"res" }
    # Check input
    comn.check_input(params, default_params, critical_params)  


    number_of_states = params["number_of_states"]
    istep   = params["isnap"]
    fstep   = params["fsnap"]
    outdir = params["outdir"]
    verbosity = params["verbosity"]


    ci_coefficients_libra = []
    nSDs = len( sd_states_unique_sorted[0] ) + 1

    # Add one to the number of CI states because the ground state is not included yet
    nCIs  = number_of_states + 1

    if verbosity >=1:
        print("=============== In make_T_matrices ============") 
        print(F"number_of_states = {number_of_states} nCIs = {nCIs} nSDs = {nSDs}")

    # Results
    SD2CI = []

    # Compute
    for step in range( fstep - istep ):
        SD2CI.append( CMATRIX( nSDs, nCIs ) )

        # Make the list of ci_coefficients for each step in the way Libra accepts
        # ci_coefficients_libra.append( [] ) -  deprecate AVA 7/8

        # Start with the ground state. This is not explicitly given by electronic strcture calculations
        # Deprecate - AVA 7/8
        #ci_coefficients_libra[step].insert( 0, [0.0] * nSDs )
        #ci_coefficients_libra[step][0][0] = 1.0 

        # Ground state "CI"
        SD2CI[step].set(0, 0, 1.0+0.0j)

        if verbosity >=2:
            print(F"   ========== in step {step} ===========")
            print(F"    sd_states_unique_sorted{step} ", sd_states_unique_sorted[step] )

        # Excited states "CI"s
        for i in range( number_of_states ):
            #count = 0

            # The ci wavefunction is a linear combination of SD states. Make a list of zeros the size of the number of unique
            # SD states + 1 for the ground state
            #ci_coefficients_libra[step].append( [0.0] * nSDs ) - deprecate - AVA 7/8


            # For each ci_coefficient in this ci state for this step, get the ci coefficients and spin (alp or bet)
            nconf = len(ci_coefficients[step][i])  

            if verbosity >= 3:
                print(F"      handling CI {i}, number of configurations = {nconf},  ci_coeffs" ,ci_coefficients[step][i] )
            

            # Exclude ground state here in the index, that info is not explicitly contained 
            # in the ci_coefficients_dynamics list from electronic structure calculations
            tmp_ci_basis_state_and_spin = []
            for k in range( nconf ):
                tmp_ci_basis_state_and_spin.append( [ci_basis_states[step][i][k] , spin_components[step][i][k]] )

            if verbosity >= 3:
                print("      tmp_ci_basis_state_and_spin = ", tmp_ci_basis_state_and_spin)


            # Now, loop over the SDs (excluding the ground state) to assign the coefficients
            for j in range( nSDs-1 ):

                # Check to see if one of the SDs from the list of unique SDs comprises this ci state
                if sd_states_unique_sorted[step][j] in tmp_ci_basis_state_and_spin:

                    # ok, it has found a match, now what is the index of the SD in the list of unique SDs?
                    item_index = tmp_ci_basis_state_and_spin.index(sd_states_unique_sorted[step][j])
                    val = float(ci_coefficients[step][i][item_index]) * (1.0+0.0j)

                    if verbosity >= 4:
                        print(F"      found match: j = {j}, what = { sd_states_unique_sorted[step][j]}, its index = {item_index} ")                        
                        print(F" {val} goes to matrix element {i+1}, {j+1}")
                    
                    #ci_coefficients_libra[step][i+1][j+1] = float(ci_coefficients[step][i][item_index])  # deprecate - AVA 7/8
                    SD2CI[step].set(j+1, i+1, float(ci_coefficients[step][i][item_index]) * 1.0+0.0j )


        if verbosity >= 2: 
            print("Matrix before the normalization")
            data_outs.print_matrix(SD2CI[step])

        # Sanity check. Make sure sum of squared elements of columns == 1:
        for i in range( nCIs ):
            check_norm = 0.0
            for j in range( nSDs ):
                #check_norm += ci_coefficients_libra[step][i][j]**2
                check_norm += abs(SD2CI[step].get(j, i))**2
            check_norm = 1.0 / math.sqrt(check_norm)

            SD2CI[step].scale(-1, i, check_norm*(1.0+0.0j))
           

        if verbosity >= 2: 
            print("Matrix after the normalization")
            data_outs.print_matrix(SD2CI[step])

            
        # Output the transformation matrix. This is how you can double check that it worked ( it does ... :) )
        SD2CI[step].show_matrix( "%s/T_%s.txt" % (outdir, str(step)) )

    return SD2CI





def make_T_matrices_fast( ci_coefficients, ci_basis_states, spin_components, sd_basis_states_unique, _params):
    """
    This function makes the "T"ransformation matricies that convert between the SD basis to the CI-like (or many-body (MB)) basis.

    This funciton is made to be used within the NBRA Libra workflow, where things such as ci_coefficients, ci_basis_states, spin_components, 
    and sd_states_unique_sorted have been extracted from TD-DFT calculations. As of 11/30/2020, compatable ES programs
    include CP2K, DFTB+ and Gaussian.

    Args:
        ci_coefficients (list of lists of lists): coefficients for the many-body states for each step
        ci_basis_states (list of lists): All SD basis states that comprise the many-body excitations for each step
        spin_components (list of lists): the spin components of the excitation (alpha or beta excitaiton?) for all states and all steps  
        sd_basis_states_unique (list): 1 of each of the SP transitions (and its spin) that made up the considered CI states
        _params (dict): control parameters
            * **_params["number_of_states"]** (int): number of excited MB states [default: 1]
            * **_params["isnap"] (int): step at which to start counting  [ default: 0 ]
            * **_params["fsnap"] (int): step at which to stop counting [ default: 1]
            * **_params["outdir"] (string): output directory for the T matricies [ default: "res"]
            * **_params["verbosity"] (int): level of extra m messages you wanna see [ default: 0 - no debug info ] 

    Returns:
        SD2CI (list of CMATRIX): CMATRIX at each timestep where the rows are SDs and the cols are MB states. The columns contain the coefficients of the MB expansion for each MB state

    """

    params = dict(_params)

    critical_params = []
    # Default parameters
    default_params = { "number_of_states":1, "isnap":0, "fsnap":1, "verbosity":0, "outdir":"res" }
    # Check input
    comn.check_input(params, default_params, critical_params)  


    number_of_states = params["number_of_states"]
    istep   = params["isnap"]
    fstep   = params["fsnap"]
    outdir = params["outdir"]
    verbosity = params["verbosity"]


    ci_coefficients_libra = []
    nSDs = len( sd_basis_states_unique ) + 1

    # Add one to the number of CI states because the ground state is not included yet
    nCIs  = number_of_states + 1

    if verbosity >=1:
        print("=============== In make_T_matrices ============") 
        print(F"number_of_states = {number_of_states} nCIs = {nCIs} nSDs = {nSDs}")
        print(F"sd_basis_states_unique ", sd_bais_states_unique )


    # Results
    SD2CI = []

    # Compute
    for step in range( fstep - istep ):

        if verbosity >=2:
            print(F"   ========== in step {step} ===========")

        SD2CI.append( CMATRIX( nSDs, nCIs ) )

        # Ground state "CI"
        SD2CI[step].set(0, 0, 1.0+0.0j)

        # Excited states "CI"s
        for i in range( number_of_states ):

            # For each ci_coefficient in this ci state for this step, get the ci coefficients and spin (alp or bet)
            nconf = len(ci_coefficients[step][i])  

            if verbosity >= 3:
                print(F"      handling CI {i}, number of configurations = {nconf},  ci_coeffs" ,ci_coefficients[step][i] )
            

            # Exclude ground state here in the index, that info is not explicitly contained 
            # in the ci_coefficients_dynamics list from electronic structure calculations
            for k in range( nconf ):
                conf = [ ci_basis_states[step][i][k] , spin_components[step][i][k] ]

                if conf in sd_basis_states_unique:    
                    indx = sd_basis_states_unique.index(conf)
                    val = ci_coefficients[step][i][k] * (1.0+0.0j)

                    SD2CI[step].set(indx+1, i+1, val )


        if verbosity >= 2: 
            print("Matrix before the normalization")
            data_outs.print_matrix(SD2CI[step])

        # Sanity check. Make sure sum of squared elements of columns == 1:
        for i in range( nCIs ):
            check_norm = 0.0
            for j in range( nSDs ):
                #check_norm += ci_coefficients_libra[step][i][j]**2
                check_norm += abs(SD2CI[step].get(j, i))**2
            check_norm = 1.0 / math.sqrt(check_norm)

            SD2CI[step].scale(-1, i, check_norm*(1.0+0.0j))
           

        if verbosity >= 2: 
            print("Matrix after the normalization")
            data_outs.print_matrix(SD2CI[step])

            
        # Output the transformation matrix. This is how you can double check that it worked ( it does ... :) )
        SD2CI[step].show_matrix( "%s/T_%s.txt" % (outdir, str(step)) )

    return SD2CI





def run(_params):
    """
    This function utilizes the single-particle results, results of TD-DFT(B) calculations, and
    user's definition of excited state basis to sompute the the Hvib in the many-body basis

    Recent update by AVA:
    if we are working with the SI states, we don't can about the ordering of SDs, so in the stage of forming T matrices, we
    can remove a lot of critical redundancies. We will be using the sorting type "identity" in that case (now default)


    Args:
       _params (dict): control parameters

    Returns:
       None
    """

    params = dict(_params)

    critical_params = []
    # Default parameters
    path = os.getcwd()
    default_params = {  "data_set_paths": [ F"{path}/res/"], "logfile_directory":F"{path}/all_logfiles",
                        "read_S_data" : 1, "read_S_re":1,  "read_S_im":0,
                        "S_data_re_prefix": "S_ks_",  "S_data_re_suffix": "_re",
                        "S_data_im_prefix": "S_ks_",  "S_data_im_suffix": "_im",
                        "read_St_data" : 1, "read_St_re":1,  "read_St_im":0,
                        "St_data_re_prefix": "St_ks_",  "St_data_re_suffix": "_re",
                        "St_data_im_prefix": "St_ks_",  "St_data_im_suffix": "_im",
                        "read_hvib_data" : 1, "read_hvib_re":1,  "read_hvib_im":0,
                        "hvib_data_re_prefix": "E_ks_",  "hvib_data_re_suffix": "_re",
                        "hvib_data_im_prefix": "E_ks_",  "hvib_data_im_suffix": "_im",
                        "isnap":0, "fsnap":1, 
                        "data_dim":1, "active_space":[0], "orbital_indices":[0], "homo_index":0,
                        "es_software":"dftb+", "isUKS":0, "tolerance":0.01, "number_of_states":1,
                        "orbital_normalization":False, "orbital_phase_correction":False,
                        "state_normalization":True, "state_phase_correction":True,
                        "do_state_reordering":2, "state_reordering_alpha":0,
                        "sorting_type":"identity", "dt": 1.0*units.fs2au,
                        "outdir":F"{path}/res_mb_sp/"
                     }

    # Check input
    comn.check_input(params, default_params, critical_params)  


    homo_index = params["homo_index"]
    orbital_indices = params["orbital_indices"]
    res_dir = params["outdir"]
    dt = params["dt"]
    orbital_normalization = params["orbital_normalization"]
    orbital_phase_correction = params["orbital_phase_correction"]
    state_normalization = params["state_normalization"]
    state_phase_correction = params["state_phase_correction"]
    start_time = params["isnap"]
    sorting_type = params["sorting_type"]


    log_file = open("step3_many_body.log", "w")
    log_file.close()

    #================== Preparation ================
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)


    #=============== Single-particle stuff ================
    # 1.1. Get the single-particle information       
    S, St, E  = step3.get_step2_data(params)    

    # 1.2. Orthogonalize orbitals and apply the phase correction to them - experimental
    ndata_sets = len(S)
    for i in range(ndata_sets):
        if( orbital_normalization ):
            step3.apply_normalization( S[i], St[i] )
        if( orbital_phase_correction ):
            step3.apply_phase_correction( St[i] )


    #=============== TD-DFT(B) stuff ================    
    # 2.1. Get the information from the TDDFT calculations
    res = get_step2_mb_sp_properties( params )

    # Unpack variables
    sd_basis_states_unique = res[0] 
    ci_basis_states = res[1]
    ci_coefficients = res[2]
    ci_energies = res[3]
    spin_components = res[4]    


    with open("step3_many_body.log", "a") as log_file:
        print("===============", file = log_file)
        print( len(sd_basis_states_unique),  file = log_file )
        print("SD basis states unique", file = log_file )
        for sd in sd_basis_states_unique:
            print ( sd , file = log_file )



    # 2.2. Reindex the single-particle excitations into a format expected by Libra: keeping in mind Slater Determinants notation
    sd_states_reindexed = step2_many_body.reindex_cp2k_sd_states( homo_index, orbital_indices,  sd_basis_states_unique, sd_format=2 )
    

    res1 = step3.sort_unique_SD_basis( E[0], sd_basis_states_unique, sd_states_reindexed, params )

    E_sd = res1[0]
    sd_states_unique_sorted = res1[1]
    sd_states_reindexed_sorted = res1[2]
    reindex_nsteps = res1[3]

    with open("step3_many_body.log", "a") as log_file:
        print("===============", file = log_file)
        print( len(sd_states_reindexed),  file = log_file )
        print("Reindexed SD states", file = log_file )
        for sd in sd_states_reindexed:
            print ( sd , file = log_file )

        print("===============", file = log_file)
        print( len(sd_states_unique_sorted),  file = log_file )
        print("Unique sorted SD states", file = log_file )
        for sd in sd_states_unique_sorted:
            print ( sd , file = log_file )

        print("===============", file = log_file)
        print( len(sd_states_reindexed_sorted),  file = log_file )
        print("Reindexed SD states", file = log_file )
        for sd in sd_states_reindexed_sorted:
            print ( sd , file = log_file )


    #=================== SD-level stuff =================
    # 3.1. Slater-determinant calculations: compute overlaps and time-overlaps between the SD states

    nsteps = params["fsnap"] - params["isnap"]

    # FIXME:  these instructions are for a single data set so far
    # 3.1.1. Overlaps
    S_sd = []
    for step in range(nsteps):
        # Re-indexed SD basis functions for bra and ket SDs
        bra = sd_states_reindexed_sorted[step]
        ket = sd_states_reindexed_sorted[step]
        s_sd = mapping.ovlp_mat_arb( bra, ket, S[0][step], False )
        S_sd.append( s_sd )
        s_sd.real().show_matrix("%s/S_sd_%d_re" % (res_dir, int(start_time+step)))

   
    # 3.1.2. Time-overlaps
    St_sd = []
    for step in range(nsteps-1):
        # Re-indexed SD basis functions for bra and ket SDs
        bra = sd_states_reindexed_sorted[step]
        ket = sd_states_reindexed_sorted[step+1]
        st_sd = mapping.ovlp_mat_arb( bra, ket, St[0][step], False )
        St_sd.append( st_sd )
        st_sd.real().show_matrix("%s/St_sd_%d_re" % (res_dir, int(start_time+step)))
    

    #==================== CI-level stuff ========================
    # 4.1. Read the the list of excitation energies at each timestep, and compute the midpoints
    ci_midpoint_energies = compute_ci_energies_midpoint( ci_energies, params )

    # 4.2. Print out the T matrices
    if sorting_type=="identity":
        # The new, faster way
        SD2CI = make_T_matrices_fast( ci_coefficients, ci_basis_states, spin_components, sd_basis_states_unique,  params )
    else:
        # The old way
        SD2CI = make_T_matrices( ci_coefficients, ci_basis_states,  spin_components, sd_states_unique_sorted,  params )

    # 4.3. Compute overlaps and time-overlaps in the CI basis
    S_ci, St_ci  = [], []
    for step in range( nsteps ):
        s_ci  = SD2CI[step].H() * S_sd[step]  * SD2CI[step]
        S_ci.append(  s_ci  )

    for step in range( nsteps-1 ):
        st_ci = SD2CI[step].H() * St_sd[step] * SD2CI[step+1]
        St_ci.append( st_ci )



    #==================== Applying corrections ======================
    # 5.1. Apply orthonormalization, state reordering, and phase correction to the many-body basis
    if( state_normalization ):
        step3.apply_orthonormalization_general( S_ci, St_ci )

    step3.apply_state_reordering_general( St_ci, ci_midpoint_energies, params )

    if( state_phase_correction ):
        step3.apply_phase_correction_general( St_ci )
 

    # 5.2. Save S and St for CI states to disk
    print("Output the CI data to the res directory..." )
    for step in range( nsteps ):
        S_ci[step].real().show_matrix("%s/S_ci_%d_re" % (res_dir, int(start_time+step)))
    for step in range( nsteps-1 ):
        St_ci[step].real().show_matrix("%s/St_ci_%d_re" % (res_dir, int(start_time+step)))

    # 5.3. Make the Hvib in the many-body basis
    Hvib = [ [] ]
    ci_hvib = None
    for step in range( nsteps-1 ): 
        ci_nacs = (  0.5j / dt ) * CMATRIX ( ( St_ci[step] - St_ci[step].H() ).real() )    
        ci_hvib = ci_midpoint_energies[step] - ci_nacs
        Hvib[0].append( ci_hvib)
        ci_hvib.real().show_matrix("%s/Hvib_ci_%d_re" % (res_dir, int( start_time+step )))
        ci_hvib.imag().show_matrix("%s/Hvib_ci_%d_im" % (res_dir, int( start_time+step )))
    Hvib[0].append( ci_hvib)    # appending the last element twice to make it nsteps


    return Hvib


