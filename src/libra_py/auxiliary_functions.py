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
from libra_py import cube_plotting_methods
from libra_py import CP2K_methods
from libra_py.workflows.nbra import mapping
from libra_py.workflows.nbra import step3
from libra_py import units

# Concurrency for large systems
import concurrent.futures








# This file is temp only. These functions will eventually be placed in Libra somewhere ...
def curr_and_final_step_job( istep, fstep, njobs, njob ):

    total_steps = fstep - istep + 1

    nsteps_per_job = int(total_steps/njobs)

    if njob == 0:
        curr_step  = istep
        final_step = istep + nsteps_per_job - 1
        print(curr_step, final_step)

    elif njob > 0 and njob < ( njobs - 1 ):
        curr_step  = istep + njob * nsteps_per_job - 1
        final_step = istep + ( njob + 1 ) * nsteps_per_job - 1
        print(curr_step, final_step)

    elif njob == ( njobs - 1 ):
        curr_step  = istep + njob * nsteps_per_job - 1
        final_step = fstep
        print(curr_step, final_step)

    return curr_step, final_step




def cp2k_distribute( istep, fstep, nsteps_this_job, cp2k_positions, cp2k_input, curr_job_number ):
    """
    Distributes cp2k jobs for trivial parallelization 

    Make sure that your cp2k input file has absolute paths to the following input parameters:
        BASIS
        POTENTIAL
        < any other file dependent parameter (dftd3, etc) >

    Args:
        
    """

    # Now we need to distribute the jobs into job batches
    # First, make a working directory where the calculations will take place
    os.chdir("wd")

    nsteps = fstep - istep + 1
    njobs  = int( nsteps / nsteps_this_job ) 

    curr_step = istep

    os.system("mkdir job"+str(curr_job_number)+"")
    os.chdir("job"+str(curr_job_number)+"")

    os.system("cp ../../"+cp2k_positions+" .")
    os.system("cp ../../"+cp2k_input+" .")

    #os.system("cp ../../submit_template.slm .")
    # Now, we need to edit the submit file

    # Now, in jobs folder njob, we should do only a certain number of steps
    for step in range( nsteps_this_job ):

        read_trajectory_xyz_file( cp2k_positions, curr_step )

        # Now, we need to edit the cp2k_input file
        tmp = open(cp2k_input)
        A   = tmp.readlines()
        sz  = len(A)
        tmp.close()
        tmp2 = open("step_%d"%curr_step+".inp","w"); tmp2.close()
        for i in range(sz):

            b = A[i].strip().split()

            if not b:
                continue

            tmp2 = open("step_%d"%curr_step+".inp","a")

            if b[0].lower() == "PROJECT".lower() or b[0].lower() == "PROJECT_NAME".lower():
                tmp2.write("  %s %s \n" % ("PROJECT", "step_%d"%curr_step +"_sp"))
            elif b[0].lower() == "COORD_FILE_NAME".lower():
                tmp2.write("  %s %s \n" % ("COORD_FILE_NAME", "../../cp2k_atomic_coordinates/coord_%d"%curr_step+".xyz"))
            else:
                tmp2.write(A[i])
            tmp2.close()

        curr_step += 1

    os.chdir("../../") 




def read_trajectory_xyz_file(file_name: str, t: int):
    """
	This function reads the trajectory of a molecular dynamics .xyz file and
	extract the 't' th step then writes it to coord-t.xyz file. This function
	is used in single point calculations for calculations of the NACs where one has
	previously obtained the trajectory via any molecular dynamics packages. 
	
	Args:
	
	    file_name (string): The trajectory .xyz file name.
		
		t (integer): The desired time to extract its .xyz 
		             coordinates which starts from zero.
	
	"""
	
    f = open(file_name,'r')
    lines = f.readlines()
    f.close()
	
    # The number of atoms for each time step in the .xyz file of the trajectory.
    number_of_atoms = int(lines[0].split()[0])

    # Write the coordinates of the 't' th step in file coord-t.xyz
    f = open('coord-%d'%t+'.xyz','w')

    # This is used to skip the first two lines for each time step.
    n = number_of_atoms+2

    for i in range( n * t, n * ( t + 1 ) ):
        f.write( lines[i] )
    f.close()





def CP2K_input_static(cp2k_sample_energy_input: str, project_name: str,\
                           trajectory_file: str, time_step: int):
    """
	This function provides the inputs required for single point calculations used for 
	calculation of NACs and overlap matrix. This function requires a sample CP2K input
	used for ENERGY calculations and the project name which will used to identify the 
	produced cube files. 
	
	Args:
	
	    cp2k_sample_energy_input (string): The CP2K sample input file for energy calculations.
		
		project_name (string): The poject name of the CP2K sample input (in &GLOBAL-> PROJECT)
		                       This is required for identifying the cube files.
		
		trajectory_file (string): The .xyz trajectory file obtained from previous 
		                          molecular dynamics calculations
								  
		time_step (integer): The time step of the single point calculation. It starts from 0.
		
	"""
    
    # Reading the .xyz trajectory file and extract the coordinates at time 'time_step'.
    read_trajectory_xyz_file(trajectory_file, time_step)
	
	# The new project name for each time_step.
    f_new = project_name+'-%d'%time_step+'.inp'
	
    f = open(cp2k_sample_energy_input,'r')
    lines = f.readlines()
    f.close()

    # wfn_restart_file name
    c = 0

	
    f = open(f_new,'w')
    for i in range(0,len(lines)):
	    # Writing the COORD_FILE_NAME in the input
        if 'COORD_FILE_NAME'.lower() in lines[i].lower().split():
            f.write('      COORD_FILE_NAME')
            f.write(' coord-%d'%time_step+'.xyz')
            f.write('\n')
		# Writing the project name
        elif 'PROJECT'.lower() in lines[i].lower().split():
            f.write('  PROJECT')
            f.write(' %s'%project_name+'-%d'%time_step)
            f.write('\n')
        elif 'PROJECT_NAME'.lower() in lines[i].lower().split():
            f.write('  PROJECT')
            f.write(' %s'%project_name+'-%d'%time_step)
            f.write('\n')
        elif 'WFN_RESTART_FILE_NAME'.lower() in lines[i].lower().split():
            f.write('    WFN_RESTART_FILE_NAME')
            f.write(' %s'%project_name+'-%d'%(time_step-1)+'-RESTART.wfn')
            f.write('\n')
        #elif 'WFN_RESTART_FILE_NAME'.lower() in lines[i].lower().split():
        #    if c==0:
        #        f.write('    WFN_RESTART_FILE_NAME')
        #        f.write(' %s'%project_name+'-%d'%(time_step-1)+'-RESTART.tdwfn')
        #        f.write('\n')
        #    elif c==1:
        #        f.write('    WFN_RESTART_FILE_NAME')
        #        f.write(' %s'%project_name+'-%d'%(time_step-1)+'-RESTART.wfn')
        #        f.write('\n')
        #    c = c + 1
        elif 'SCF_GUESS'.lower() in lines[i].lower().split():
            f.write('    SCF_GUESS RESTART')
            f.write('\n')
        else:
            f.write(lines[i])




def normalize_ci_coefficients(ci_coefficients_raw_unnorm):
    """
    This funciton normalizez a list of lists of numbers
    """
    
    nstates = len(ci_coefficients_raw_unnorm)
    ci_coefficients_raw_norm = []
       
    for i in range(nstates):
        
        norm = 0.0
        ci_coefficients_raw_norm.append( [] )
        
        for j in ci_coefficients_raw_unnorm[i]:
            
            norm += j*j
        
        norm = math.sqrt(norm)
        
        # numpy way
        ci_coefficients_raw_unnorm[i] = np.array(ci_coefficients_raw_unnorm[i])
        ci_coefficients_raw_norm[i] = abs ( ci_coefficients_raw_unnorm[i] / norm )
        ci_coefficients_raw_norm[i] = list(ci_coefficients_raw_norm[i])
        
        # non-numpy, pythonic way
        #ci_coefficients_raw_norm[i] = [ ci_coefficients_raw_unnorm[i][k] / norm for k in range(len(ci_coefficients_raw_unnorm[i])) ]
    
    return ci_coefficients_raw_norm
        



def get_es_output(params):
    """
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
    """

    overlap_matrix = np.zeros( ( len( cubefiles_prev ) , len( cubefiles_curr ) ) )

    for i in range( 0, len( cubefiles_prev ) ):
        for j in range( 0, len( cubefiles_curr ) ):
            overlap_matrix[i][j] = cube_file_methods.integrate_cube( cubefiles_prev[i], cubefiles_curr[j], dv )

    return overlap_matrix



def compute_cube_ks_overlaps( cubefiles_prev, params):
    """
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




def plot_cubes( params ):
    """
    Plotting the cubes for selected energy levels using VMD
    """

    # Unpack the Kohn-Sham orbital indicies and the Kohn-Sham orbitals to be plotted. Also unpack the 
    # path to the directory where the molecular orbitals will be plotted
    ks_orbital_indicies  = params["ks_orbital_indicies"]
    states_to_be_plotted = params["states_to_be_plotted"]
    # For VMD
    path_to_tcl_file = params["path_to_tcl_file"]
    # The molecular orbital images directory
    MO_images_directory  = params["MO_images_directory"]

    # isUKS flag for spin-polarized and spin-unpolarized
    isUKS  = int( params["isUKS"] )
    # The current step
    curr_step = int(params["curr_step"])
    # If the path does not exist create it.
    if not os.path.isdir(MO_images_directory):
        os.makedirs(MO_images_directory)

    # Extracting the phase factor from params
    phase_factor_visual = params["phase_factor_visual"]

    min_band = ks_orbital_indicies[0]

    cubefile_names_prev = CP2K_methods.cube_file_names_cp2k( params )

    tcl_file = open(path_to_tcl_file,'r')
    tcl_lines = tcl_file.readlines()
    tcl_file.close()


    if isUKS == 1:
        alp_cubefile_names_prev = cubefile_names_prev[0::2]
        phase_factor_alpha      = phase_factor_visual[0::2]

        bet_cubefile_names_prev = cubefile_names_prev[1::2]
        phase_factor_beta       = phase_factor_visual[1::2]

        for state_to_be_plotted in states_to_be_plotted:
            alpha_cube_name = alp_cubefile_names_prev[state_to_be_plotted-min_band]
            beta_cube_name = bet_cubefile_names_prev[state_to_be_plotted-min_band]

            new_file_alpha = open("vmd_alpha_cube_plot_%d.tcl" % curr_step,'w')
            new_file_beta = open("vmd_beta_cube_plot_%d.tcl" % curr_step,'w')

            for j in range(len(tcl_lines)):
                if 'mol load cube' in tcl_lines[j]:
                    new_file_alpha.write( 'mol load cube %s\n' % alpha_cube_name )
                    new_file_beta.write( 'mol load cube %s\n' % beta_cube_name )
                elif 'render TachyonInternal' in tcl_lines[j]:
                    new_file_alpha.write( 'render TachyonInternal %s/%s.tga\n' % ( MO_images_directory, alpha_cube_name.replace('cubefiles/','').replace('.cube','') ) )
                    new_file_beta.write( 'render TachyonInternal %s/%s.tga\n' % ( MO_images_directory, beta_cube_name.replace('cubefiles/','').replace('.cube','') ) )
                elif 'isosurface' in tcl_lines[j].lower():
                    tmp_elements_alpha = tcl_lines[j].split()
                    tmp_elements_alpha[5] = str( phase_factor_alpha[state_to_be_plotted-min_band] * float( tmp_elements_alpha[5] ) )
                    isosurface_line_alpha = ' '.join(tmp_elements_alpha)
                    new_file_alpha.write( isosurface_line_alpha + '\n' )


                    tmp_elements_beta = tcl_lines[j].split()
                    tmp_elements_beta[5] = str( phase_factor_beta[state_to_be_plotted-min_band] * float( tmp_elements_beta[5] ) )
                    isosurface_line_beta = ' '.join(tmp_elements_beta)
                    new_file_beta.write( isosurface_line_beta + '\n' )

                else:
                    new_file_alpha.write( tcl_lines[j] )
                    new_file_beta.write( tcl_lines[j] )

            new_file_alpha.close()
            new_file_beta.close()

            os.system('vmd < vmd_alpha_cube_plot_%d.tcl' % curr_step)
            os.system('vmd < vmd_beta_cube_plot_%d.tcl' % curr_step)
            #os.system('rm vmd_alpha_cube_plot_%d.tcl' % curr_step)
            #os.system('rm vmd_beta_cube_plot_%d.tcl' % curr_step)

    else:

        for state_to_be_plotted in states_to_be_plotted:
            cube_name = cubefile_names_prev[state_to_be_plotted-min_band]

            new_file = open("vmd_cube_plot_%d.tcl" % curr_step,'w')
            for j in range(len(tcl_lines)):
                if 'mol load cube' in tcl_lines[j]:
                    new_file.write( 'mol load cube %s\n' % cube_name )
                elif 'render TachyonInternal' in tcl_lines[j]:
                    new_file.write( 'render TachyonInternal %s/%s.tga\n' % ( MO_images_directory,cube_name.replace('cubefiles/','').replace('.cube','') )  )
                elif 'isosurface' in tcl_lines[j].lower():
                    tmp_elements = tcl_lines[j].split()
                    tmp_elements[5] = str( phase_factor_visual[state_to_be_plotted-min_band] * float( tmp_elements[5] ) )
                    isosurface_line = ' '.join(tmp_elements)
                    new_file.write( isosurface_line + '\n' )
                else:
                    new_file.write( tcl_lines[j] )

            new_file.close()

            os.system('vmd < vmd_cube_plot_%d.tcl' % curr_step)
            #os.system('rm vmd_cube_plot_%d.tcl' % curr_step)





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



def read_energies_from_cp2k_log_file( cp2k_log_file_name: str, time: int, min_band: int, max_band: int, spin: int ):
    """
    This function reads the energies from CP2K log file for a specific spin.
    
    Args:
    
        cp2k_log_file_name (str): The output file produced by CP2K
        
        time (int): The time step of the molecular dynamics. For single point calculations
                    it should be set to 0.
                    
        min_band (int): The minimum state to be considered.
        
        max_band (int): The maximum state to be considered.
        
        spin (int): This number can only accept two values: 1 for alpha and 2 for beta spin.
        
    Returns:
    
        E (1D numpy array): The vector consisting of the KS energies from min_band to max_band.
        
    """
    # First openning the file and reading its lines
    f = open( cp2k_log_file_name, 'r' )
    lines = f.readlines()
    f.close()
    
    # The lines containing 'Eigenvalues of the occupied subspace'
    lines_with_occupied   = []
    # The lines containing 'Eigenvalues of the unoccupied subspace'
    lines_with_unoccupied = []
    
    # The lines containing the energies of the occupied states
    occ_energies_last_line   = []
    # The lines containing the energies of the unoccupied states
    unocc_energies_last_line = []
    
    total_energy = 0.0

    for i in range(0,len(lines)):
        # Find the occupied states lines
        if 'Eigenvalues of the occupied subspace'.lower() in lines[i].lower():
        
            if  str( spin ) in lines[i].lower().split():
                lines_with_occupied.append( i )
                for j in range( i, len( lines ) ):
                    # Find the last line number containing the energies for occupied states
                    if len( lines[j].split() ) == 0:
                        occ_energies_last_line.append( j )
                        break
        # Find the unoccupied states lines
        if 'Eigenvalues of the unoccupied subspace'.lower() in lines[i].lower():
            if  str( spin ) in lines[i].lower().split():
                lines_with_unoccupied.append(i)
                for j in range( i, len( lines ) ):
                    # Find the last line number containing the energies for unoccupied states
                    if len( lines[j].split() ) == 0:
                        unocc_energies_last_line.append( j )
                        break
       
        if "Total energy:" in lines[i]:  
            total_energy = float( lines[i].split()[2] )


    print("Debugging ....")
    print( lines_with_occupied[time]      )
    print( occ_energies_last_line[time]   )
    print( lines_with_unoccupied[time]    )
    print( unocc_energies_last_line[time] )

    # The Kohn-Sham energies
    ks_energies = []
    # Start after two lines of the 'Eigenvalues of the occupied subspace'
    for i in range( lines_with_occupied[time] + 2, occ_energies_last_line[time] - 1 ):
        for j in range(0,len(lines[i].split())):
            ks_energies.append(float(lines[i].split()[j]))
    
    # Start after two lines of the 'Eigenvalues of the unoccupied subspace'
    for i in range( lines_with_unoccupied[time] + 2, unocc_energies_last_line[time] ):
        if not 'Reached'.lower() in lines[i].lower().split():
            for j in range(0,len(lines[i].split())):
                ks_energies.append(float(lines[i].split()[j]))

    # Now appending all the energies into a numpy array
    ks_energies = np.array(ks_energies)
        
    # Returning the energeis from min_band to max_band
    return ks_energies[min_band-1:max_band], total_energy




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




def run_cp2k_job_number( params ):
    """
    """

    ks_orbital_indicies = [ int(i) for i in list( params["ks_orbital_indicies"].split(",") ) ]
    params.update({"ks_orbital_indicies":ks_orbital_indicies})
    ks_orbital_homo_index = int(params["ks_orbital_homo_index"]) 
    nsteps_this_job = int(params["nsteps_this_job"])
    cp2k_exe = params["cp2k_exe"]
    cp2k_input_template = params["cp2k_input_template"]
    project_name = params["project_name"]
    trajectory_xyz_filename = params["trajectory_xyz_filename"]
    njob    = int(params["njob"])
    nprocs  = int(params["nprocs"])
    res_dir = params["res_dir"]
    dt = float(params["dt"])

    logfile_directory = params["logfile_directory"]
    es_software = params["es_software"]
    number_ci_states = int(params["number_ci_states"])
    tolerance = float(params["tolerance"])    

    output_level = int(params["output_level"])

    do_phase_corrections = int(params["do_phase_corrections"])

    curr_step = int(params["istep"])
    params.update({ "curr_step":curr_step })

    do_cube_visualization = int(params["do_cube_visualization"])


    # For empty lists for the KS, SD, and CI bases. These are to be used for collecting data at each step
    S_ks_job  = []
    St_ks_job = []
    E_ks_job  = []
    total_energies_job = [] 
 
    S_sd_job  = []
    St_sd_job = []
    sd_basis_states_unique = []

    S_ci_job  = []
    St_ci_job = []
    ci_coefficients_job   = []
    ci_basis_states_job   = []
    ci_energies_job       = []
    spin_components_job   = []

    # Initial step

    # Set up the cp2k input template with the atomic positions for the first timestep of this job batch
    CP2K_input_static( cp2k_input_template, project_name, trajectory_xyz_filename, curr_step )
    timer1 = time.time()
    os.system("mpirun -np %d %s -i %s-%i.inp -o logfiles/step_%d.log "%( nprocs, cp2k_exe, project_name, curr_step, curr_step ) )
    os.system("mv *.pdos pdosfiles")
    os.system("mv *.cube cubefiles")
    print("CP2K time was ", time.time() - timer1)

    # We have compted the first SCF calculation for this job, now to read the output data and cubes
    print("\nReading the initial step using pool")
    pool = mp.Pool( processes = nprocs )
    print("\nJust called pool")
    cubefile_names_prev = CP2K_methods.cube_file_names_cp2k( params )
    print("\nInitial step for this job, prinring cubefile_names_prev\n", cubefile_names_prev)   
    #sys.exit(0)
    cubefiles_prev = pool.map( cube_file_methods.read_cube, cubefile_names_prev )
    pool.close()

    # If we are to visualize the cubefiles, then we should compute the grid-mesh before hand to save time later
    if do_cube_visualization == 1:
        phase_factor_visual = np.ones( ( ks_orbital_indicies[-1] - ks_orbital_indicies[0] + 1 ) * 2 )
        params.update({"phase_factor_visual":phase_factor_visual})
        states_to_be_plotted = [ int(i) for i in list( params["states_to_be_plotted"].split(",") ) ]
        params.update({"states_to_be_plotted":states_to_be_plotted})
        plot_cubes( params )

    # After reading the cubefiles for curr_step, we should delete them to be memory efficient
    os.system("rm cubefiles/*")

    hvib_ks_re, total_energy = form_Hvib_real( "logfiles/step_"+str(params["curr_step"])+".log", 0, ks_orbital_indicies[0], ks_orbital_indicies[-1], int(params["isUKS"]))
    E_ks_re = data_conv.nparray2CMATRIX( hvib_ks_re )
    E_ks_job.append( E_ks_re )
    total_energies_job.append( total_energy )

    # Get the es output - this was originally coded for the CI basis. So we need another function for the ks only
    excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = get_es_output( params )
    ci_coefficients_raw_norm = normalize_ci_coefficients(ci_coefficients_raw_unnorm)

    # Extract the uniquie SD basis states from the ci basis states
    for ci_basis_state_index in range( len( ci_basis_raw ) ):
        for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):
           sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ]
           if sd_basis_state_and_spin not in sd_basis_states_unique:
               sd_basis_states_unique.append( sd_basis_state_and_spin )
    print( "\nSD basis states = ", sd_basis_states_unique )

    ci_basis_states_job.append( ci_basis_raw )
    ci_coefficients_job.append( ci_coefficients_raw_norm )
    ci_energies_job.append( excitation_energies )
    spin_components_job.append( spin_components )
    print("\n###############")
    print("Printing cp2k TD-DFPT logfile information for step number", curr_step)
    print("Excitation energies\n", excitation_energies)
    print("CI basis raw\n", ci_basis_raw)
    print("Spin components\n", spin_components)
    print("Unnormalized CI coefficients\n", ci_coefficients_raw_unnorm)
    print("Normalized CI coefficients\n", ci_coefficients_raw_norm)

    curr_step += 1

    # All other steps after initial step for this job
    for step in range( nsteps_this_job-1 ):

        params.update({ "curr_step":curr_step })
        CP2K_input_static( cp2k_input_template, project_name, trajectory_xyz_filename, curr_step )
        timer1 = time.time()
        os.system("mpirun -np %d %s -i %s-%i.inp -o logfiles/step_%d.log "%( nprocs, cp2k_exe, project_name, curr_step, curr_step ) )
        os.system("mv *.pdos pdosfiles")
        os.system("mv *.cube cubefiles")
        print("CP2K time was ", time.time() - timer1)

        ########################################################
        ##### Collect stuff from log and cube files here  ######
        ########################################################

        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = get_es_output( params )
        ci_coefficients_raw_norm = normalize_ci_coefficients(ci_coefficients_raw_unnorm)
        print("\n###############")
        print("Printing cp2k TD-DFPT logfile information for step number", curr_step)
        print("Excitation energies\n", excitation_energies)
        print("CI basis raw\n", ci_basis_raw)
        print("Spin components\n", spin_components)
        print("Unnormalized CI coefficients\n", ci_coefficients_raw_unnorm)
        print("Normalized CI coefficients\n", ci_coefficients_raw_norm)

        hvib_ks_re, total_energy = form_Hvib_real( "logfiles/step_"+str(params["curr_step"])+".log", 0, ks_orbital_indicies[0], ks_orbital_indicies[-1], int(params["isUKS"]))
        E_ks_re = data_conv.nparray2CMATRIX( hvib_ks_re )
        E_ks_job.append( E_ks_re )
        total_energies_job.append( total_energy )

        # Now, read in the cube files for steps nstep-1 and nstep
        # and form the S_ks and St_ks overlap matricies

        cubefiles_curr, S_ks_prev, S_ks_curr, St_ks = compute_cube_ks_overlaps( cubefiles_prev, params )
        cubefiles_prev = cubefiles_curr

        if do_cube_visualization == 1:
            for row_index in range(len(St_ks)):
                if St_ks[row_index][row_index]<0:
                    phase_factor_visual[row_index] = phase_factor_visual[row_index] * (-1)
            params.update({"phase_factor_visual":phase_factor_visual})
            plot_cubes( params )
        #sys.exit(0)

        os.system("rm cubefiles/*")
        #sys.exit(0)

        # Now, we need to orthonormalize the KS basis. Recall that in the function step3.normalization, 
        # only the St variables are updated, not the S variables. 
        S_ks_prev = data_conv.nparray2CMATRIX( S_ks_prev )
        S_ks_curr = data_conv.nparray2CMATRIX( S_ks_curr )
        St_ks     = data_conv.nparray2CMATRIX( St_ks )
        step3.apply_normalization( [S_ks_prev, S_ks_curr] , [St_ks] )

        S_ks_job.append(S_ks_prev)
        St_ks_job.append(St_ks)
        if step == nsteps_this_job-2:
            S_ks_job.append(S_ks_curr)

        # Extract the uniquie SD basis states from the ci basis states
        for ci_basis_state_index in range( len( ci_basis_raw ) ):
            for sd_basis_state_index in range( len( ci_basis_raw[ ci_basis_state_index ] ) ):
               sd_basis_state_and_spin = [ ci_basis_raw[ ci_basis_state_index ][ sd_basis_state_index ] , spin_components[ ci_basis_state_index ][ sd_basis_state_index ] ] 
               if sd_basis_state_and_spin not in sd_basis_states_unique:
                   sd_basis_states_unique.append( sd_basis_state_and_spin )
        print( "\nSD basis states = ", sd_basis_states_unique )

        ci_basis_states_job.append( ci_basis_raw )
        ci_coefficients_job.append(   ci_coefficients_raw_norm )
        ci_energies_job.append( excitation_energies )
        spin_components_job.append( spin_components )

        curr_step += 1
        print("\nJust finished with step", step)      
 
    print("\nWe just finished with all of the step")

    # Now we have the S_ks, St_ks and E_ks for each step in this job folder. We now output the Hvib on the KS level
    # Output KS data to res directory
    print("\nLength of S_ks_job and St_ks_job:", len(S_ks_job), len(St_ks_job))
    for step in range( nsteps_this_job ):
        S_ks_job[step].real().show_matrix("%s/S_ks_%d_re" % (res_dir, int(params["istep"])+step))
        E_ks_job[step].real().show_matrix("%s/E_ks_%d_re" % (res_dir, int(params["istep"])+step))

    for step in range( nsteps_this_job-1 ):
        St_ks_job[step].real().show_matrix("%s/St_ks_%d_re" % (res_dir, int(params["istep"])+step))
    #sys.exit(0)

    # Now, time to compute S_sd and St_sd
    # Start by reindexing the unique Slater determinant basis. The current SD bases are not able to be read by Libra
    sd_states_reindexed = reindex_cpk2_sd_states( ks_orbital_homo_index, ks_orbital_indicies, sd_basis_states_unique )
    print("\nsd_state_reindex = \n", sd_states_reindexed)
    #sys.exit(0)

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

        print("\nWe have to test to see if the reindexing is working")
        print("\n The energy of the SD for this job is") 
        print(e)
        print("\n The reindex for this job is")  
        print(reindex)
        print("\n The old sd list")         
        print(sd_states_reindexed)
        print("\n The new sd list")         
        print(sd_states_reindexed_sorted[step])
        print("\n sd_states_unique_sorted")
        print(sd_states_unique_sorted[step])
    #sys.exit(0)

    # For each step make S_sd and St_sd
    for step in range( nsteps_this_job ):
        s_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step], S_ks_job[step],  use_minimal=False )
        S_sd_job.append( s_sd )
    for step in range( nsteps_this_job-1 ):
        st_sd = step3.mapping.ovlp_mat_arb(  sd_states_reindexed_sorted[step], sd_states_reindexed_sorted[step+1], St_ks_job[step],  use_minimal=False )
        St_sd_job.append( st_sd )
    #sys.exit(0)

    # Output Slater determinant data to the res directory
    for step in range( nsteps_this_job ):
        S_sd_job[step].real().show_matrix("%s/S_sd_%d_re" % (res_dir, int(params["istep"])+step))
        E_sd_job[step].real().show_matrix("%s/E_sd_%d_re" % (res_dir, int(params["istep"])+step))
    for step in range( nsteps_this_job-1 ):
        St_sd_job[step].real().show_matrix("%s/St_sd_%d_re" % (res_dir, int(params["istep"])+step))
    #sys.exit(0)

    # Now, we have computed the overlaps in the KS and SD bases. Now, we need to compute the SD2CI matrix for each step
    # Now, form the SD2CI matrix for each step
    ci_coefficients = []
    # Add one to the number of CI states because the ground state is not included yet
    nSDs = len( sd_basis_states_unique ) + 1
    nCIs = number_ci_states + 1
    SD2CI = []
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

            tmp_ci_basis_state_and_spin = [ [ci_basis_states_job[step][i][k] , spin_components_job[step][i][k]] for k in range(len(ci_coefficients_job[step][i])) ]
            print("\n step ", step, "ci state ", i+1)
            print("\n sd_basis_states_in_this_ci_state ", tmp_ci_basis_state_and_spin) 
            print("\n sd_states_unique_sorted[step]    ", sd_states_unique_sorted[step])
            print("\n ci_coefficients_job[step]        ", ci_coefficients_job[step])
            #sys.exit(0)

            for j in range( nSDs-1 ):
                if sd_states_unique_sorted[step][j] in tmp_ci_basis_state_and_spin:   
                    # ok, it has found a match, now what is the index?
                    item_index = tmp_ci_basis_state_and_spin.index(sd_states_unique_sorted[step][j])
                    print("\nitem_index", item_index)
                    ci_coefficients[step][i+1][j+1] = float(ci_coefficients_job[step][i][item_index])    
        print("\n ci_coefficients[step] ", ci_coefficients[step])
 
        SD2CI.append( CMATRIX( nSDs, nCIs ) )
        for i in range( nSDs ):
            for j in range( nCIs ):
                SD2CI[step].set( i, j, ci_coefficients[step][j][i] * (1.0+0.0j) )
        SD2CI[step].show_matrix("T_"+str(step)+".txt")
    #sys.exit(0)

    # For each step make S_ci and St_ci
    for step in range( nsteps_this_job ):
        s_ci = SD2CI[step].H()  * S_sd_job[step]  * SD2CI[step]
        S_ci_job.append( s_ci )
    for step in range( nsteps_this_job-1 ):
        st_ci = SD2CI[step].H() * St_sd_job[step] * SD2CI[step+1]
        St_ci_job.append( st_ci )

    # Now, compute the CI energy matrix at each-point and the mid-points
    # For each step
    ci_energies_job_cmatrix = []
    for step in range( nsteps_this_job ):
        ci_energies_job_cmatrix.append( CMATRIX( number_ci_states + 1, number_ci_states + 1 ) )
        for state in range( number_ci_states + 1 ):
            if state == 0:
                ci_energies_job_cmatrix[step].set( state, state, total_energies_job[step] )
            else:
                ci_energies_job_cmatrix[step].set( state, state, total_energies_job[step] + ( ci_energies_job[step][state-1]  * units.ev2Ha )  )

    # At the midpoints
    print("\n Excitation energies = \n", ci_energies_job)
    ci_midpoint_energies = []
    for step in range( nsteps_this_job-1 ):
        total_energy_mid_point = 0.5 * ( total_energies_job[step] + total_energies_job[step+1] )
        ci_midpoint_energies.append( CMATRIX( number_ci_states + 1, number_ci_states + 1 ) )
        for state in range( number_ci_states + 1 ):
            if state == 0:
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point )
            else:
                midpoint_energy = 0.5 * ( ci_energies_job[step][state-1] + ci_energies_job[step+1][state-1] )
                ci_midpoint_energies[step].set( state, state, total_energy_mid_point + ( midpoint_energy  * units.ev2Ha )  )
    #sys.exit(0)

    # Are we to perform state reordering?
    if int(params["perform_state_reordering"]) == 1:
        params2 = {"do_state_reordering":int(params["do_state_reordering"]), "state_reordering_alpha":float(params["state_reordering_alpha"])}
        print("\nApplying state reordering")
        step3.apply_state_reordering_ci( St_ci_job, ci_midpoint_energies, params2 )
    #sys.exit(0)

    # Are we to perform phase corrections?
    if do_phase_corrections == 1:
        print("\nApplying phase corrections")
        step3.apply_phase_correction_ci( St_ci_job )
    #sys.exit(0)

    # Output CI data to res directory
    print("\nLength of S_ci_job and St_ci_job:", len(S_ci_job), len(St_ci_job))
    for step in range( nsteps_this_job ):
        S_ci_job[step].real().show_matrix("%s/S_ci_%d_re"   % (res_dir, int(params["istep"])+step))
        ci_energies_job_cmatrix[step].real().show_matrix("%s/E_ci_%d_re"   % (res_dir, int(params["istep"])+step))
    for step in range( nsteps_this_job-1 ):
        St_ci_job[step].real().show_matrix("%s/St_ci_%d_re" % (res_dir, int(params["istep"])+step))

    # Now, compute the CI NACs and compute the CI Hvib
    for step in range( nsteps_this_job-1 ): 
        ci_nacs = -(  0.5j / dt ) * CMATRIX ( ( St_ci_job[step] - St_ci_job[step].H() ).real() )    
        ci_hvib = ci_midpoint_energies[step] + ci_nacs
        ci_hvib.real().show_matrix("%s/Hvib_ci_%d_re" % (res_dir, int(params["istep"])+step))
        ci_hvib.imag().show_matrix("%s/Hvib_ci_%d_im" % (res_dir, int(params["istep"])+step))

