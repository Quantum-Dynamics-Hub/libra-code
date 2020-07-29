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
.. module:: CP2K_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing the CP2K inputs and outputs.
.. moduleauthors:: 
       Mohammad Shakiba, Brendan Smith, Alexey V. Akimov 
  
"""


import os
import sys
import math
import re
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *




def ndigits( integer_number: int ):
    """
    This function calculates the number of digits for an integer number.

    Args:

        integer_number ( integer ): An integer number.

    Returns:
 
        digit_num ( integer ): The number of digits for 'num'.
    """

    # Setting up the digit number counter
    digit_num = 0

    while( integer_number > 0 ):

        digit_num = digit_num+1
        integer_number = integer_number//10
        
    return digit_num



def state_num_cp2k(state_num: int):
    """
    This function turns an integer number to a string. This string 
    has the format in five digits and if the number of digits is less than
    five it will append zeros before the number. This is necessary to read
    the .cube files produced by CP2K.

    Args:

        state_num (integer): Integer number which is in fact the energy level.

    Returns:

        state_num_str (string): The number in five digits (or more) in string format.

    """
    # Setting up an initial string
    state_num_str = ''

    # Now count the number of digits for the state sumber
    # and then add 0 before the number in the string

    if ndigits( state_num ) == 1:

        state_num_str = '0000%d' % state_num

    elif ndigits( state_num ) == 2:

        state_num_str = '000%d' % state_num

    elif ndigits( state_num ) == 3:

        state_num_str = '00%d' % state_num

    elif ndigits( state_num ) == 4:

        state_num_str = '0%d' % state_num

    else:
        state_num_str = '%d' % state_num
    
    return state_num_str




def cube_file_names_cp2k( static_or_dynamics, project_name, min_band, max_band, time_step, spin ):
    """
    This function will produce the cube file names produced by CP2K both for energy calculations
    and molecular dynamics. 

    Args:

        static_or_dynamics ( integer ): This variable has two options:
                                      0: static calculations
                                      1: molecular dynamics calculations

        project_name ( string ): The project name used in CP2K input file.

        min_band ( integer ): The minimum state number to be considered.

        max_band ( integer ): The maximum state number to be considered.

        time_step ( integer ): The time step in CP2K calculations.

        spin ( integer ): This variable can get only two options:
                          1 for alpha spin
                          2 for beta spin

    Returns:

        names ( list ): The list of cube file names.
	
    """
    # Define the names of the cube files for time t.
    cube_file_names = []

    # Dynamics calculations
    if static_or_dynamics == 1:

        for i in range( min_band, max_band + 1 ):

            # Using the state_num_cp2k function to obtain the names of the .cube files for state 'i'.
            state_num = state_num_cp2k( i )

            # For time step '( t )'
            cube_file_name_of_band_i = str( '%s-WFN_%s_%d-1_%d.cube' % ( project_name, state_num, spin, time_step ) )
            cube_file_names.append( cube_file_name_of_band_i )

    # Static calculations
    elif static_or_dynamics == 0:

        for i in range( min_band, max_band+1 ):

            # Using the int_to_five_digit_str to obtain the names of the .cube files for state 'i'.
            state_num = state_num_cp2k( i )

            # For time step '( t )'
            cube_file_name_of_band_i = str('%s-%d-WFN_%s_%d-1_0.cube' % ( project_name, time_step, state_num, spin ) )
            cube_file_names.append( cube_file_name_of_band_i )
			
    return cube_file_names



def read_cp2k_tddfpt_log_file( tddfpt_log_file_name, num_ci_states, tolerance ):
    """
    This function reads log files generated from TD-DFPT calculations using CP2K and returns the TD-DFPT
    excitation energies, Slater determinent states in terms of the Kohn-Sham orbital indicies, 
    and the Slater determinent coefficients that comprise the multi-configurational electronic states
    
    Args:

        tddfpt_log_file_name ( string ): name of the log file for this particular timestep

        num_ci_states ( int ): how many ci states to consider

        tolerance ( float ): the tolerance for accepting SDs are members of the CI wavefunctions

    Returns:

       excitation_energies ( list ): The excitation energies in the log file.

       ci_basis ( list ): The CI-basis configuration.

       ci_coefficients ( list ): The coefficients of the CI-states.

    """
    
    f = open(tddfpt_log_file_name,'r')
    lines = f.readlines()
    f.close()

    for i in range( 0, len(lines) ):

        if 'excitation' in lines[i].lower().split():
            if 'analysis' in lines[i].lower().split():

                # When found the line in which contains 'Excitation analysis' in 
                # the log file, append it in the variable exc_anal_line
                exc_anal_line = i
        
        if 'states' in lines[i].lower().split():
            if 'multiplicity' in lines[i].lower().split():

                # Here we search for the line that contains
                # 'R-TDDFPT states of multiplicity 1' in the log file
                r_tddfpt_line = i

    excitation_energies = []
    
    # Start from 5 lines after finding the line contaning 'R-TDDFPT states of multiplicity 1'
    # This is because they contain the energies from that line.
    for i in range( r_tddfpt_line+5, len( lines ) ):
        if len( lines[i].split() ) == 0:
            break
        excitation_energies.append( float( lines[i].split()[2] ) )
        
    # Start from 5 lines after finding the line contaning 'Excitation analysis'
    # From that point we have the state numbers with their configurations.
    # So, we append the lines which contain only 'State number' and stop
    # whenever we reach to a blank line.
    state_num_lines = []
    for i in range( exc_anal_line+5, len( lines ) ):

        if len( lines[i].split() ) == 0:
            state_num_lines.append( i )
            break

        if len( lines[i].split() ) == 1:
            state_num_lines.append( i )

    # Setting up the CI-basis list and their coefficients list
    ci_basis        = []
    ci_coefficients = []
    for i in range( num_ci_states ):
        
        # CI states and their coefficients for each excited state
        tmp_ci_state              = []
        tmp_ci_state_coefficients = []
        for j in range( state_num_lines[i]+1, state_num_lines[i+1] ):

            if abs( float( lines[j].split()[2] ) ) > tolerance:
                tmp_ci_state.append( [ int( lines[j].split()[0] ), int( lines[j].split()[1] ) ]  )
                tmp_ci_state_coefficients.append( float( lines[j].split()[2] ) )

        # Append the CI-basis and and their coefficients for
        # this state into the ci_basis and ci_coefficients lists
        ci_basis.append( tmp_ci_state )
        ci_coefficients.append( tmp_ci_state_coefficients )
        
    return excitation_energies, ci_basis, ci_coefficients



def cp2k_distribute( istep, fstep, nsteps_this_job, cp2k_positions, cp2k_input, curr_job_number ):
    """
    Distributes cp2k jobs for trivial parallelization 

    *** Make sure that your cp2k input file has absolute paths to the following input parameters:
        BASIS
        POTENTIAL
        < any other file dependent parameter (dftd3, etc) >

    Args:
    
        istep (int): The initial time step for molecular dynamics trajectory which 
	             the user wants to start the calculations with. This paramter 
		     is set in the input file for submitting the jobs.
        
	fstep (int): The final time step for molecular dynamics trajectory which 
	             the user wants to start the calculations with. This paramter 
		     is set in the input file for submitting the jobs.
		     
        nsteps_this_job (int): The number of step for the job to be distributed.
	
	cp2k_positions (str): The full path to the molecular dynamics trajectory xyz file.
	
	curr_job_number (int): The current job number.
        
    """

    # Now we need to distribute the jobs into job batches
    # First, make a working directory where the calculations will take place
    os.chdir("wd")

    # total number of steps
    nsteps = fstep - istep + 1
    # number of jobs obtained from nsteps/nsteps_this_job
    njobs  = int( nsteps / nsteps_this_job ) 
    # The current step is set to istep
    curr_step = istep
    # Making the job directories in the wd folder like job0, job1, ...
    os.system("mkdir job"+str(curr_job_number)+"")
    # Chnage directory to job folder
    os.chdir("job"+str(curr_job_number)+"")
    # Copy the trajectory xyz file
    os.system("cp ../../"+cp2k_positions+" .")
    # Copy the CP2K sample input
    os.system("cp ../../"+cp2k_input+" .")

    #os.system("cp ../../submit_template.slm .")
    # Now, we need to edit the submit file

    # Now, in job folder njob, we should do only a certain number of steps
    for step in range( nsteps_this_job ):
        # extract the curr_step xyz coordianates from the trajectory file and write it to another xyz file
        read_trajectory_xyz_file( cp2k_positions, curr_step )

        # Now, we need to edit the cp2k_input file
        tmp = open(cp2k_input)
        A   = tmp.readlines()
        sz  = len(A)
        tmp.close()
	# create a new input file based on the curr step
        tmp2 = open("step_%d"%curr_step+".inp","w"); tmp2.close()
        for i in range(sz):

            b = A[i].strip().split()

            if not b:
                continue

            tmp2 = open("step_%d"%curr_step+".inp","a")
            # Write the project name in the new input file
            if b[0].lower() == "PROJECT".lower() or b[0].lower() == "PROJECT_NAME".lower():
                tmp2.write("  %s %s \n" % ("PROJECT", "step_%d"%curr_step +"_sp"))
            # Wite the coordinate file name obtained for the curr_step in the new input file
            elif b[0].lower() == "COORD_FILE_NAME".lower():
                tmp2.write("  %s %s \n" % ("COORD_FILE_NAME", "../../cp2k_atomic_coordinates/coord_%d"%curr_step+".xyz"))
            else:
                tmp2.write(A[i])
            tmp2.close()

        curr_step += 1
    # Go back to the main directory
    os.chdir("../../") 


	
def read_trajectory_xyz_file(file_name: str, step: int):
    """
	This function reads the trajectory of a molecular dynamics .xyz file and
	extract the 'step' th step then writes it to coord-step.xyz file. This function
	is used in single point calculations for calculations of the NACs where one has
	previously obtained the trajectory via any molecular dynamics packages. 
	
	Args:
	
	    file_name (string): The trajectory .xyz file name.
		
            step (integer): The desired time to extract its .xyz 
		         coordinates which starts from zero.
	
    """
	
    f = open(file_name,'r')
    lines = f.readlines()
    f.close()
	
    # The number of atoms for each time step in the .xyz file of the trajectory.
    number_of_atoms = int(lines[0].split()[0])

    # Write the coordinates of the 't' th step in file coord-t.xyz
    f = open('coord-%d'%step+'.xyz','w')

    # This is used to skip the first two lines for each time step.
    n = number_of_atoms+2

    # Write the coordinates of the 'step'th time step into the file
    for i in range( n * step, n * ( step + 1 ) ):
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
        # Writing the WFN_RESTART_FILE_NAME of the previous time step
	elif 'WFN_RESTART_FILE_NAME'.lower() in lines[i].lower().split():
            f.write('    WFN_RESTART_FILE_NAME')
            f.write(' %s'%project_name+'-%d'%(time_step-1)+'-RESTART.wfn')
            f.write('\n')
        # Setting the scf guess to restart, CP2K will automatically switch to atomic if the wfn file does not exist
        elif 'SCF_GUESS'.lower() in lines[i].lower().split():
            f.write('    SCF_GUESS RESTART')
            f.write('\n')
        else:
            f.write(lines[i])



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
    # Set the total energy to zero
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




