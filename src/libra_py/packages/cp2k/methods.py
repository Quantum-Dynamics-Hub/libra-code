#*********************************************************************************
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
import time

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn   

from libra_py import data_outs
from libra_py import units



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




def cube_file_names_cp2k( params ):
    """
    This function will produce the cube file names produced by CP2K both for energy calculations
    and molecular dynamics. 

    Args:

        params( dictionary ): contains dynamical paramters
	
            project_name ( string ): The project name used in CP2K input file.
	    
            min_band ( integer ): The minimum state number to be considered.
	    
            max_band ( integer ): The maximum state number to be considered.
	    
            curr_step ( integer ): The time current step in CP2K calculations.
	    
            isUKS ( int ): 1 if spin-polarized calculations used

    Returns:

        names ( list ): The list of cube file names.
	
    """
    # Critical parameters
    critical_params = [ "project_name", "min_band", "max_band", "curr_step" ]
    # Default parameters
    default_params = { "isUKS":0 }
    # Check input
    comn.check_input(params, default_params, critical_params) 
  
    project_name = params["project_name"]
    min_band  = params["min_band"]
    max_band  = params["max_band"]
    curr_step = params["curr_step"]
    isUKS = params["isUKS"]
    es_software = params["es_software"] 

    # isUKS needs to be epressed as a string, even though it is used as a Boolean in other parts of the workflows
    spin = []
    if isUKS == 1:
        spin = [1,2]
    else:
        spin = [1]

    # Define the names of the cube files for time t.
    cube_file_names = []   

    for i in range( min_band, max_band+1 ):

        # Using the function state_num_cp2k to obtain the names of the .cube files for state 'i'.
        state_num = state_num_cp2k( i )
 
        # For the current time step "curr_step"
        if isUKS == 1:

            if es_software == "cp2k" or es_software == "gaussian":

                cube_file_name_of_band_i = str('cubefiles/%s-%d-WFN_%s_%d-1_0.cube' % ( project_name, curr_step, state_num, spin[0] ) )
                cube_file_name_of_band_j = str('cubefiles/%s-%d-WFN_%s_%d-1_0.cube' % ( project_name, curr_step, state_num, spin[1] ) )
                cube_file_names.append( cube_file_name_of_band_i )
                cube_file_names.append( cube_file_name_of_band_j )

            elif software == "dftb+":
                print("\nReading cubefiles from spin-polarized calculations not currently available")
                print("Exiting now")
                sys.exit(0)
 
        else:

            # Using the function state_num_cp2k to obtain the names of the .cube files for state 'i'.
            state_num = state_num_cp2k( i )
            cube_file_name_of_band_i = str('cubefiles/%s-%d-WFN_%s_%d-1_0.cube' % ( project_name, curr_step, state_num, spin[0] ) )
            cube_file_names.append( cube_file_name_of_band_i )
		
    return cube_file_names





def read_cp2k_tddfpt_log_file( params ):
    """
    This function reads log files generated from TD-DFPT calculations using CP2K and returns the TD-DFPT
    excitation energies, Slater determinent states in terms of the Kohn-Sham orbital indicies, 
    and the Slater determinent coefficients that comprise the multi-configurational electronic states
    
    Args:
        params ( dictionary ): parameters dictionary 
	
            logfile_name ( string ): name of the log file for this particular timestep
	    
            number_of_states ( int ): how many ci states to consider
	    
            tolerance ( float ): the tolerance for accepting SDs are members of the CI wavefunctions
	    
            isUKS ( Boolean ): flag for whether or not a spin-polarized Kohn-Sham basis is being used. TRUE means
                                 that a spin-polarized Kohn-Sham basis is being used.
    Returns:
        excitation_energies ( list ): The excitation energies in the log file.
	
        ci_basis ( list ): The CI-basis configuration.
	
        ci_coefficients ( list ): The coefficients of the CI-states.
	
        spin_components (list): The spin component of the excited states.
	
    """

    # Critical parameters
    critical_params = [ "logfile_name", "number_of_states" ]
    # Default parameters
    default_params = { "tolerance":0.05, "isUKS": 0}
    # Check input
    comn.check_input(params, default_params, critical_params)
    
    logfile_name = params["logfile_name"]
    number_of_states = int(params["number_of_states"])
    tolerance = float(params["tolerance"])
    isUKS = int(params["isUKS"])


    f = open( logfile_name, 'r' )
    lines = f.readlines()
    f.close()

    for i in range( 0, len(lines) ):
        tmp_line = lines[i].lower().split()
        if 'excitation' in tmp_line:
            if 'analysis' in tmp_line:

                # When found the line in which contains 'Excitation analysis' in 
                # the log file, append it in the variable exc_anal_line
                exc_anal_line = i

        if 'states' in tmp_line:
            if 'multiplicity' in tmp_line:

                # Here we search for the line that contains 'R-TDDFPT states of multiplicity 1'
                # or 'U-TDDFPT states of multiplicity 1' in the log file
                r_tddfpt_line = i

    excitation_energies = []
    # Old version
    ## Start from 5 lines after finding the line contaning 'R-TDDFPT states of multiplicity 1'
    ## This is because they contain the energies from that line.
    ##for i in range( r_tddfpt_line+5, len( lines ) ):
    ##    tmp_line = lines[i].split()
    ##    if len( tmp_line ) == 0:
    ##        break
    ##    excitation_energies.append( float( tmp_line[2] ) )

    for i in range(len(lines)):
        if 'TDDFPT|' in lines[i]:
            try:
                tmp = lines[i].split()
                exc_ener = float(tmp[2])
                excitation_energies.append(exc_ener)
            except:
                pass

    # Start from 5 lines after finding the line contaning 'Excitation analysis'
    # From that point we have the state numbers with their configurations.
    # So, we append the lines which contain only 'State number' and stop
    # whenever we reach to a blank line.
    state_num_lines = []
    for i in range( exc_anal_line+5, len( lines ) ):
        tmp_line = lines[i].split()
        if isUKS == 1:

            if len( tmp_line ) == 0 or '----' in lines[i]:
                state_num_lines.append( i )
                break
            if len( tmp_line )==1 or len( tmp_line )==3:
                state_num_lines.append(i)

        else:

            if len( tmp_line ) == 0 or '----' in lines[i]:
                state_num_lines.append( i )
                break
            if len(lines[i].split())==1:
                state_num_lines.append(i)
            elif len( tmp_line ) > 1:
                if tmp_line[0].isdigit() and not( tmp_line[1].isdigit() ):
                    state_num_lines.append( i )


    # Setting up the CI-basis list and their coefficients list
    ci_basis = []
    ci_coefficients = []
    spin_components = []
    for i in range( number_of_states ):

        # states and their coefficients for each excited state
        tmp_state              = []
        tmp_state_coefficients = []
        tmp_spin = []
        for j in range( state_num_lines[i]+1, state_num_lines[i+1] ):
            
            # Splitting lines[j] into tmp_splitted_line
            tmp_splitted_line = lines[j].split()

            # If we are using the spin-polarized Kohn-Sham basis, then the size of the split line is different
            # from the size of the split line when using the spin-unpolarized Kohn-Sham basis
            if isUKS == 1:

                ci_coefficient = float( tmp_splitted_line[4] )
                if ci_coefficient**2 > tolerance:
                    # We need to remove the paranthesis from the 2nd element of the temporary splitted line
                    tmp_spin.append( tmp_splitted_line[1].replace('(','').replace(')','') )
                    tmp_state.append( [ int( tmp_splitted_line[0] ), int( tmp_splitted_line[2] ) ]  )
                    tmp_state_coefficients.append( ci_coefficient  )

            # Here, we have the spin-unpolarize Kohn-Sham basis
            # For this case, spin-components will just return all alpha
            else:
                ci_coefficient = float( tmp_splitted_line[2] )
                if ci_coefficient**2 > tolerance:
                    tmp_spin.append( "alp" )
                    tmp_state.append( [ int( tmp_splitted_line[0] ), int( tmp_splitted_line[1] ) ]  )
                    tmp_state_coefficients.append( ci_coefficient  )

        # Append the CI-basis and and their coefficients for
        # this state into the ci_basis and ci_coefficients lists
        ci_basis.append( tmp_state )
        ci_coefficients.append( tmp_state_coefficients )
        spin_components.append( tmp_spin )

    return excitation_energies[0:number_of_states], ci_basis, ci_coefficients, spin_components




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
    
    Returns:
    
        None
	
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
    
    Returns:
    
        None
	
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





def CP2K_input_static(cp2k_sample_energy_input: str, project_name: str,
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
    
    Returns:
    
        None
	
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
        elif 'WFN_RESTART_FILE_NAME'.lower() in lines[i].lower().split() and ".wfn" in lines[i].lower() and "!" not in lines[i].lower():	
            f.write('    WFN_RESTART_FILE_NAME')	
            f.write(' %s'%project_name+'-%d'%(time_step-1)+'-RESTART.wfn')	
            f.write('\n')	
        # Writing the WFN_RESTART_FILE_NAME of the previous time step for the tddft calculation	
        elif 'WFN_RESTART_FILE_NAME'.lower() in lines[i].lower().split() and ".tdwfn" in lines[i].lower() and "!" not in lines[i].lower():	
            f.write('    WFN_RESTART_FILE_NAME')	
            f.write(' %s'%project_name+'-%d'%(time_step-1)+'-RESTART.tdwfn')
            f.write('\n')
        # Setting the scf guess to restart, CP2K will automatically switch to atomic if the wfn file does not exist
        elif 'SCF_GUESS'.lower() in lines[i].lower().split():
            f.write('    SCF_GUESS RESTART')
            f.write('\n')
        else:
            f.write(lines[i])



def read_energies_from_cp2k_log_file( params ):
    """
    This function reads the energies from CP2K log file for a specific spin.
    
    Args:
    
        params (dictionary):
    
            logfile_name (str): The output file produced by CP2K
        
            time (int): The time step of the molecular dynamics. For single point calculations
                        it should be set to 0.
                    
            min_band (int): The minimum state number.

            max_band (int): The maximum state number.
 
            spin (int): This number can only accept two values: 1 for alpha and 2 for beta spin.
        
    Returns:
    
        E (1D numpy array): The vector consisting of the KS energies from min_band to max_band.
        
        total_energy (float): The total energy obtained from the log file.
        
    """
    
    # Critical parameters
    critical_params = [ "logfile_name", "min_band", "max_band" ]
    # Default parameters
    default_params = { "time":0, "spin": 1}
    # Check input
    comn.check_input(params, default_params, critical_params)
    
    cp2k_log_file_name = params["logfile_name"]
    # Time step, for molecular dynamics it will read the energies
    # of time step 'time', but for static calculations the time is set to 0
    time = params["time"]
    
    # Koh-Sham orbital indicies
    # ks_orbital_indicies = params["ks_orbital_indicies"]
    
    # The minimum state number
    min_band = params["min_band"] # ks_orbital_indicies[0]
    # The maximum state number
    max_band = params["max_band"] # ks_orbital_indicies[-1]
    
    spin = params["spin"]
    
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
  

def read_energies_from_cp2k_md_log_file( params ):
    """
    This function reads the energies from CP2K molecular dynamics log file for a specific spin. The difference between this function
    and the function `read_energies_from_cp2k_log_file` is that it is more efficient and faster for molecular dynamics energy levels.
    
    Args:
    
        params (dictionary):
    
            logfile_name (str): The output file produced by CP2K
        
            init_time (int): The initial time step of the molecular dynamics. 
                        
            final_time (int): The final time step of the molecular dynamics. 
                    
            min_band (int): The minimum state number.

            max_band (int): The maximum state number.
 
            spin (int): This number can only accept two values: 1 for alpha and 2 for beta spin.
        
    Returns:
    
        E (1D numpy array): The vector consisting of the KS energies from min_band to max_band.
        
        total_energy (float): The total energy obtained from the log file.
        
    """
    
    # Critical parameters
    critical_params = [ "logfile_name", "min_band", "max_band" ]
    # Default parameters
    default_params = { "spin": 1, "init_time": 0, "final_time": 1}
    # Check input
    comn.check_input(params, default_params, critical_params)
    
    cp2k_log_file_name = params["logfile_name"]
    # Time step, for molecular dynamics it will read the energies
    # of time step 'time', but for static calculations the time is set to 0
    init_time = params["init_time"]
    final_time = params["final_time"]
    
    # Koh-Sham orbital indicies
    # ks_orbital_indicies = params["ks_orbital_indicies"]
    
    # The minimum state number
    min_band = params["min_band"] # ks_orbital_indicies[0]
    # The maximum state number
    max_band = params["max_band"] # ks_orbital_indicies[-1]
    
    spin = params["spin"]
    
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

    # All KS energies
    KS_energies = []
    for time in range(init_time,final_time):
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
        KS_energies.append(ks_energies[min_band-1:max_band])
        
    # Returning the energeis from min_band to max_band
    return KS_energies, total_energy


def read_homo_index(filename: str):
    """
    This function extract the HOMO index from CP2K log file. This index starts from 1.
    Args:

        filename (string): The full path to the log file.

    Returns:

        homo_index (integer): The HOMO index (the number of occupied orbitals).
    """
    file = open(filename,'r')
    lines = file.readlines()
    file.close()

    for i in range(len(lines)):
        if 'occupied' in lines[i].lower():
            if 'number' in lines[i].lower():
                homo_index = int(lines[i].split()[-1])

    return homo_index


def read_molog_file(filename: str):
    """
    This function reads the coefficiets of the molecular orbitals printed out
    during the MD or a single-point calculation for a structure using CP2K.
    The format of the MO coefficients in CP2K is in column format so the coeffiecients
    are written in a column. You can check this for an MOLog file printed out by CP2K.
    The number of columns in MOLog file is based on the DFT%PRINT%MO%NDIGITS. The larger the 
    number of digits (NDIGITS), the lower number of columns. This function is written in a way
    that will automatically extract all the coefficients for all the eigenvectors and returns them
    and their energies. It will also designed for reading the coefficients for MOs in different
    K-points if used by user in the CP2K input (DFT%KPOINT).
    
    Args:
    
        filename (string): The name of the MOLog ile.
        
    Returns:
    
        mo_coeffs (list): The list containing the MO coefficients for each K-point.
        
        mo_energies (liest): The list containing the energies of MOs for each K-point.
    
    """
    # first is there any k-point or not??
    file = open(filename,'r')
    lines = file.readlines()
    file.close()
    
    is_kpoint = False
    # The K-point will show itself at the beginning
    # so since the files are large we only search the first 
    # couple of lines, say first 10
    for i in range(len(lines)):
        if 'K-POINT' in lines[i]:
            is_kpoint = True
            break
    
    # set a timer
    timer = time.time()
    if is_kpoint:
        # If K-point
        print('Found K-point, will proceed reading the coefficients',
        'for each K-POINT...')
        # Lines with K-POINT
        kpoint_lines = []
        # Lines with Fermi energy
        fermi_lines = []
        for i in range(len(lines)):
            if 'K-POINT'.lower() in lines[i].lower():
                kpoint_lines.append(i)
            if 'Fermi'.lower() in lines[i].lower():
                fermi_lines.append(i)
        # All the coefficients
        mo_coeffs = []
        # All the energies
        mo_energies = []
        # Loop over all the K-points
        for i in range(len(kpoint_lines)):
            # start lines
            startl = kpoint_lines[i]
            if i==len(kpoint_lines)-1:
                # end line
                endl = len(lines)
            else:
                endl = kpoint_lines[i+1]
            # All the eigenvectors for each K-point
            eigenvectors = []
            # Lines with length less than or equal to 4
            leq_4_lines = []
            # Energies for a K-point
            energies = []
            # Find the leq_4_lines indicies
            for j in range(startl, endl):
                # lines with less than or equal to 4
                tmp = lines[j].split()
                if 0 < len(tmp)<=4:
                    if 'Fermi' not in tmp and 'gap' not in tmp:
                        leq_4_lines.append(j)
            # The lines with MO number
            mo_num_lines = leq_4_lines[0::3]
            # The lines wit energies
            energy_lines = leq_4_lines[1::3]
            # Find the energies for a K-point
            for j in energy_lines:
                tmp = lines[j].split()
                for k1 in range(4):
                    try:
                        energy = float(tmp[k1])
                        energies.append(energy)
                    except:
                        pass
            # Append the energies for this K-point
            mo_energies.append(np.array(energies))
            # Occupation lines, Needed to define the start and end
            # lines to find the coefficients
            occ_lines = leq_4_lines[2::3]
            # Now reading the coefficients
            for j in range(len(mo_num_lines)):
                start_l = occ_lines[j]+1
                if j==len(mo_num_lines)-1:
                    end_l = fermi_lines[i]
                else:
                    end_l = mo_num_lines[j+1]
                # This part with try and except make the code 
                # to consider for different number of columns
                # So the user does not need to specify the number
                # of columns or the number of atomic orbitals etc...
                for k1 in range(4,8):
                    # Each eigenvector for this K-point
                    eigenvector = []
                    for k2 in range(start_l, end_l):
                        tmp = lines[k2].split()
                        try:
                            eig_val = float(tmp[k1])
                            eigenvector.append(eig_val)
                        except:
                            pass
                    if len(eigenvector)!=0:
                        # Now this variable contains all the eigenvectors 
                        # for this K-point
                        eigenvectors.append(np.array(eigenvector))
            
            print('Done reading coefficients for K-point %d'%(i+1))
            # The mo_coeffs contains the eigenvectors. After this
            # we start for another K-point.
            mo_coeffs.append(eigenvectors)
        
    else:
        # The same procedure as above for K-point with this 
        # difference that there is no K-point
        fermi_lines = []
        for i in range(len(lines)):
            if 'Fermi'.lower() in lines[i].lower():
                fermi_lines.append(i)
        
        mo_coeffs = []
        mo_energies = []
        # start lines
        startl = 3 
        endl = fermi_lines[0] 
        eigenvectors = []
        leq_4_lines = []
        energies = []
        for j in range(startl, endl):
            # lines with less than or equal to 4
            tmp = lines[j].split()
            if 0 < len(tmp)<=4:
                if 'Fermi' not in tmp and 'gap' not in tmp:
                    leq_4_lines.append(j)

        mo_num_lines = leq_4_lines[0::3]
        energy_lines = leq_4_lines[1::3]
        for j in energy_lines:
            tmp = lines[j].split()
            for k1 in range(4):
                try:
                    energy = float(tmp[k1])
                    energies.append(energy)
                except:
                    pass
        mo_energies.append(np.array(energies))
        occ_lines = leq_4_lines[2::3]
        for j in range(len(mo_num_lines)):
            start_l = occ_lines[j]+1
            if j==len(mo_num_lines)-1:
                end_l = fermi_lines[0]
            else:
                end_l = mo_num_lines[j+1]
            for k1 in range(4,8):
                eigenvector = []
                for k2 in range(start_l, end_l):
                    tmp = lines[k2].split()
                    try:
                        eig_val = float(tmp[k1])
                        eigenvector.append(eig_val)
                    except:
                        pass
                if len(eigenvector)!=0:
                    eigenvectors.append(np.array(eigenvector))

        print('Done reading coefficients for the MOLog file')
        mo_coeffs.append(eigenvectors)
    print('Elapsed time:', time.time()-timer)
            
    return mo_energies, mo_coeffs


def read_ao_matrices(filename):
    """
    This function reads the files printed by CP2K for atomic orbital matrices data.
    If the AO_MATRICES keyword is activated in CP2K input, then all the specified 
    matrices in that section will be appened in one file. This function reads all the
    matrices and data.
    Args:
        filename (string): The name of the file that contains AO matrices data.
    Returs:
        data (list): A list that contains the AO matrices with the same order as in the file.
    """    
    
    file = open(filename,'r')
    lines = file.readlines()
    file.close()
    
    mat_lines = []
    data_lines = []
    for i in range(len(lines)):
        tmp = lines[i].split()
        if 'MATRIX' in lines[i]:
            mat_lines.append(i)
            #print(lines[i])
        if 0<len(tmp)<=4 and 'MATRIX' not in lines[i]:
            data_lines.append(i)
    data_lines = np.array(data_lines)
    data = []
    for i in range(len(mat_lines)):
        print(lines[mat_lines[i]])
        data_0 = []
        start_line = mat_lines[i]
        if i==len(mat_lines)-1:
            end_line = len(lines)-1
        else:
            end_line = mat_lines[i+1]
        #print(start_line, end_line)
        ind = np.where(np.logical_and(start_line<data_lines, data_lines<end_line))# and data_lines.all()<end_line)
        #print(ind[0])
        #print(data_lines[ind])
        for j in range(len(ind[0])):
            start_line_1 = int(data_lines[int(ind[0][j])])            
            if j==len(ind[0])-1:
                end_line_1 = int(end_line)
            else:
                end_line_1 = int(data_lines[int(ind[0][j+1])])
            #print(lines[start_line_1])
            #print(lines[end_line_1])
            for k1 in range(4):
                tmp = []
                try:
                    for k2 in range(start_line_1+1,end_line_1):
                        #print(lines[k2])
                        if len(lines[k2].split())>0:
                            #print(lines[k2])
                            tmp.append(float(lines[k2].split()[4+k1]))

                    data_0.append(tmp)
                except:
                    pass
        data.append(np.array(data_0))
        
    return data


def extract_coordinates(trajectory_xyz_file_name: str, time_step: int):
    """
    This function reads the trajectory xyz file and extract the coordinates of a
    time step.
    
    Args:
    
        trajcetory_xyz_file_name (string): The trajectory xyz file name.
        
        time_step (integer): The time step.
        
    Returns:
    
        coordinates (list): The list of the x, y, and z coordinates.
    """
    # Reading the file and its lines
    file = open(trajectory_xyz_file_name,'r')
    lines = file.readlines()
    file.close()
    # number of atoms obtianed from the first line of 
    # the xyz file
    natoms = int(lines[0].split()[0])
    # the coordinates list
    coordinates = []
    # start lines
    start_l = time_step*(natoms+2)
    # end line
    end_l = (time_step+1)*(natoms+2)
    for i in range(start_l+2,end_l):
        tmp = lines[i].split()
        name = tmp[0]
        # turn Angstrom into Bohr unit
        x = float(tmp[1])*units.Angst
        y = float(tmp[2])*units.Angst
        z = float(tmp[3])*units.Angst
        coordinates.append([name,x,y,z])
        
    return coordinates


def find_basis_set(basis_set_files_path: list, unique_atoms: list, basis_set_names: list):
    """
    This function searches the atoms basis sets for a set of unique atoms in 
    a set of basis sets files and reads them in a form so that we can create 
    a shell for atomic orbital overlap computation using libint.
    
    Args:
    
        basis_set_files_path (list): A list containing the full path to the 
                                     basis sets files (like BASIS_MOLOPT).
        unique_atoms (list): The list of unique atoms.
        
        basis_set_names (list): A list containing the name of the basis sets used
                                for each of the unique atoms. The length of this
                                list should be the same as unique_atoms list.
    
    Returns:
    
        basis_set_data (list): This lists contains the information of all the basis
                               sets. They include, the angular momentum values, 
                               the exponents, and the contraction coefficients.
    """
    # The variables for appending the l_values,
    # exponents, and contraction coefficients
    l_vals = []
    exponents = []
    coeffs = []
    basis_set_data = []
    # loop over all the basis set files
    for i in range(len(basis_set_files_path)):
        file = open(basis_set_files_path[i],'r')
        lines = file.readlines()
        file.close()
        # find the basis set for each atom in unique_atoms
        for j in range(len(unique_atoms)):
            # initializing the empty lists for each atom
            exp_atom = []
            coeff_atom = []
            l_vals_atom = []
            data_atom = []
            for k in range(len(lines)):
                tmp = lines[k].lower().split()
                # we use a try and except so that if it couldn't find the name it doesn't crash
                try:
                    if (unique_atoms[j].lower()) == tmp[0] and basis_set_names[j].lower() == tmp[1]:
                        print('Found basis set %s for atom %s in %s'%(basis_set_names[j],\
                                                                      unique_atoms[j],basis_set_files_path[i]))
                        print('Reading the data for atom %s in the basis file'%unique_atoms[j])
                        # start reading the basis set
                        # number of basis set
                        num_basis = int(lines[k+1])
                        # Lines with information about the basis set
                        # These lines are like this
                        # 3 0 1 5 2 2
                        # their order is like this:
                        # principal quantum number, minimum angular momentum value
                        # maximum angular momentum value, number of exponents
                        # the number of contractions for each of the angular momentum values
                        info_lines = []
                        info_lines.append(k+2)
                        if num_basis > 1:
                            c = 0
                            while len(info_lines) < num_basis:
                                tmp_1 = lines[k+2+c].split()
                                n_exponent = int(tmp_1[3])
                                c += n_exponent+1
                                info_lines.append(k+2+c)
                        for p in info_lines:
                            exp = []
                            coeffs = []
                            l_min = int(lines[p].split()[1])
                            l_max = int(lines[p].split()[2])
                            n_l = l_max-l_min+1
                            l_vals = []
                            c = 1
                            l_start = l_min
                            if l_min==l_max:
                                l_end = l_min+1
                            else:
                                l_end = l_max+1
                            for l_val in range(l_start, l_end):
                                l_val_rep = int(lines[p].split()[3+c])
                                c += 1
                                for pp in range(l_val_rep):
                                    l_vals.append(l_val)
                            n_exponent = int(lines[p].split()[3])
                            start_l = p+1
                            end_l = p+n_exponent+1
                            for p1 in range(start_l, end_l):
                                exp.append(float(lines[p1].split()[0]))
                            c = 0
                            for p2 in range(len(l_vals)):
                                c += 1
                                coeffs = []
                                for p1 in range(start_l, end_l):
                                    coeffs.append(float(lines[p1].split()[c]))
                                data_atom.append([l_vals[p2],exp,coeffs])
                            
                except:
                    pass
            if len(data_atom)!=0:
                # The basis_set data appends all the data for this atom
                basis_set_data.append([unique_atoms[j],data_atom])
                        
    return basis_set_data

def make_shell(coordinates: list, basis_set_data: list, is_spherical: bool):
    """
    This function makes a libint shell from the coordinates, basis_set_data from
    the find_basis_set function. The sphecrical or cartesian coordinates flags define
    the basis it needs to work in.
    """
    # setting up a counter for initializing the shell using 
    # liblibra_core.initialize_shell function
    c = 0
    # for each of the atoms
    for i in range(len(coordinates)):
        atom_name = coordinates[i][0]
        # making the coordinates into a VECTOR type to be able to use in C++
        coords_init = [coordinates[i][1], coordinates[i][2], coordinates[i][3]]
        a = VECTOR(coords_init[0], coords_init[1], coords_init[2])
        # for each of the atoms types basis set data
        for j in range(len(basis_set_data)):
            # check if the atom type is for that 
            # specific atom in the coordinates
            if atom_name==basis_set_data[j][0]:
                for k in range(len(basis_set_data[j][1])):
                    # l_value
                    l_val = basis_set_data[j][1][k][0]
                    # exponents
                    exp = Py2Cpp_double(list(basis_set_data[j][1][k][1]))
                    # contraction coefficients
                    coeff = Py2Cpp_double(list(basis_set_data[j][1][k][2]))
                    if c==0:
                        # if it is the first atom and
                        # the first basis set initialize the libint shell
                        shell = initialize_shell(l_val, is_spherical, exp, coeff, a)
                    else:
                        # for all other atoms
                        # add to the created shell 
                        add_to_shell(shell, l_val, is_spherical, exp, coeff, a)
                    # add 1 to the counter 
                    c += 1
                    
    return shell


def resort_molog_eigenvectors(l_vals):
    """
    This function returns the resotring indices for resoting the MOLog 
    eigenvectors according to this order:
    
    MOLog order (example for Cd atom):
    2s, 3s | 3py, 3pz, 3px | 4py, 4pz, 4px | 4d-2, 4d-1, 4d0, 4d+1, 4d+2 | 5d-2, 5d-1, 5d0, 5d+1, 5d+2 | ...
    However, the atomic orbital overla computed from the libint version of Psi4 is not ordered 
    as above. The ordering is like this:
    2s, 3s | 3pz, 3px, 3py | 4pz, 4px, 4py | 4d0, 4d+1, 4d-1, 4d+2, 4d-2 | 5d0, 5d+1, 5d-1, 5d+2, 5d-2 | ...
    Therefore, we need to resort the eigenvectors to be able to use the code properly.
    
    Args:
    
        l_vals (list): A list containing the angular momentum values for atoms
                       in the order of the MOLog files.
                       
    Returns:
    
        new_indices (numpy array): The new indices that needs to be used for reordering.
    
    """
    # new indices
    new_indices = []
    # setting up the counter
    c = 0
    # loop over all the angular momentum values
    for i in range(len(l_vals)):
        l_val = l_vals[i]
        # find the reorder indices needed for this l_val
        reordered_indices = index_reorder(l_val)
        for j in range(len(reordered_indices)):
            # now append it by plus the counter since
            # we aim to do it for all the eigenvectors
            # and l values
            new_indices.append(c+reordered_indices[j])
        # increase the counter with t
        c += len(reordered_indices)
    # Return the new indices
    return new_indices
    
    
def index_reorder(l_val):
    """
    This function returns the new ordering of the angular momentum value based on the 
    order used in Psi4. Detailed explanation was given for the resort_molog_eigenvectors function.
    
    Args:
    
        l_val (integer): The angular momentum value.
                       
    Returns:
    
        new_order (numpy array): The new order of the indices for the l_val.
    
    """
    
    # for s orbital
    if l_val == 0:
        new_order = [1]
    # for p orbital
    elif l_val == 1:
        new_order = [2,3,1]
    # for d orbital
    elif l_val == 2:
        new_order = [3,4,2,5,1]
    # for f orbital
    elif l_val == 3:
        new_order = [4,5,3,6,2,7,1]
    # for g orbital
    elif l_val == 4:
        new_order = [5,6,4,7,3,8,2,9,1]

    # The indeices
    return np.array(new_order)-1


def generate_translational_vectors(origin, N, periodicity_type):
    """
    This function generates the translational vectors for periodic systems.
    For monolayers the generated vectors does not add the orthogonal axis but
    for bulk all directions are added.

    Args:
        origin (list): The translational vectors are obtained with respect to this origin.
        N (list): An array that contains the number of cells to be considered in 
                     each of the X, Y, and Z directions.
        periodicity_type (string): The periodicity type. It can only get these values:
                                   'XY', 'XZ', and 'YZ' for monolayers and 'XYZ' for bulk systems.

    Returns:
        translational_vectors (numpy array): The translational vectors for that system.

    """
    # Basis vectors
    i_vec = np.array([1,0,0])
    j_vec = np.array([0,1,0])
    k_vec = np.array([0,0,1])
    # origin
    origin = np.array(origin)
    translational_vectors = []
    x_range = [0]
    y_range = [0]
    z_range = [0]
    Nx = N[0]
    Ny = N[1]
    Nz = N[2]
    if 'x' in periodicity_type.lower():
        x_range = range(-Nx,Nx+1)
    if 'y' in periodicity_type.lower():
        y_range = range(-Ny,Ny+1)
    if 'z' in periodicity_type.lower():
        z_range = range(-Nz,Nz+1)

    for n_i in x_range:
        for n_j in y_range:
            for n_k in z_range:
                vec = origin + n_i*i_vec+n_j*j_vec+n_k*k_vec
                if not np.array_equal(vec,origin):
                    translational_vectors.append(vec)
                
    return np.array(translational_vectors)


def molog_lvals(filename:str):
    """
    This function returns all the angular momentum values in the order 
    the eigenvectors are written. Unlike the molden files we need to extract 
    this from the molog files.
    
    Args:
        
        filename (string): The MOLog file name.
        
    Returns:
     
        l_vals_all (list): The list for all angular momentum values.
    """
    # Opening the file and reading all the lines
    file = open(filename,'r')
    lines = file.readlines()
    file.close()
    # all l_values including the 
    # ones which are repeated as well.
    l_vals = []
    # principal quantum number
    # this will be used to distinguish
    # between the l_values we want to read.
    p_nums = []
    # Lines with length less than or equal 
    # to one which are the molecular orbital 
    # number, their energies, and their occupation numbers
    leq_4_lines = []
    for i in range(3,len(lines)):
        tmp = lines[i].split()
        if len(tmp)<4 and len(tmp)!=0 and 'Fermi' not in lines[i] and 'MO' not in lines[i]:
            leq_4_lines.append(i)
    # We only need for one part in the MOLog file
    start_l = leq_4_lines[2]
    end_l = leq_4_lines[3]
    for i in range(start_l,end_l):
        tmp = lines[i].split()
        try:
            if 's' in tmp[3]:
                # n is the principal quantum number
                n = int(tmp[3][0])
                p_nums.append(n)
                l_vals.append(0)
            if 'p' in tmp[3]:
                n = int(tmp[3][0])
                p_nums.append(n)
                l_vals.append(1)
            if 'd' in tmp[3]:
                n = int(tmp[3][0])
                p_nums.append(n)
                l_vals.append(2)
            if 'f' in tmp[3]:
                n = int(tmp[3][0])
                p_nums.append(n)
                l_vals.append(3)
            if 'g' in tmp[3]:
                n = int(tmp[3][0])
                p_nums.append(n)
                l_vals.append(4)
        except:
            pass
    # initializing the l_vals_all
    l_vals_all = []
    # Now we need to get rid of the repeated l_values
    # Below will do it all!
    for i in range(len(l_vals)):
        if i==1:
            l_vals_all.append(l_vals[0])
        if i!=1 and i!=(len(l_vals)-1):
            if p_nums[i]!=p_nums[i-1]:
                l_vals_all.append(l_vals[i])
            if p_nums[i]==p_nums[i-1] and l_vals[i]!=l_vals[i-1]:
                l_vals_all.append(l_vals[i])
    
    return l_vals_all


def cp2k_xtb_ot_inp(ot_inp_temp: str, traj_xyz_filename: str, step: int):
    """
    This function gets the xTB orbital transformation (OT) input template and the trajectrory xyz file name
    and makes a new input based on the time step. The aim of this input is to obtain a converged
    wfn RESTART file and use it for diagonalization so that we can print out the MOLog or molden files.

    Args:

        ot_inp_temp (string): The name of the OT input template.

        traj_xyz_filename (string): The name of the trajectory xyz file name.

        step (integer): The time step.

    Returns: 

        None
    """

    file = open(ot_inp_temp)
    lines = file.readlines()
    file.close()

    read_trajectory_xyz_file(traj_xyz_filename, step)

    file = open('xtb_ot_step_%d.inp'%step,'w')

    for i in range(len(lines)):
        if 'PROJECT_NAME'.lower() in lines[i].lower():
            file.write('PROJECT_NAME OT_%d\n'%step)
        elif 'COORD_FILE_NAME'.lower() in lines[i].lower():
            file.write('COORD_FILE_NAME coord-%d.xyz\n'%step)
        elif 'WFN_RESTART_FILE_NAME'.lower() in lines[i].lower():
            file.write('WFN_RESTART_FILE_NAME OT_%d-RESTART.wfn\n'%(step-1))
        elif 'SCF_GUESS'.lower() in lines[i].lower():
            file.write('SCF_GUESS RESTART\n')
        else:
            file.write(lines[i])

    file.close()



def cp2k_xtb_diag_inp(diag_inp_temp: str, step: int):
    """
    This function gets the xTB diagonalization input template 
    and makes a new input based on the time step. This input will use the 
    wfn RESTART file name obtained from the OT calclulations.

    Args:

        diag_inp_temp (string): The name of the diagonalization input template.

        step (integer): The time step.

    Returns: 

        None
    """

    file = open(diag_inp_temp)
    lines = file.readlines()
    file.close()

    file = open('xtb_diag_step_%d.inp'%step,'w')

    for i in range(len(lines)):
        if 'PROJECT_NAME'.lower() in lines[i].lower():
            file.write('PROJECT_NAME Diag_%d\n'%step)
        elif 'COORD_FILE_NAME'.lower() in lines[i].lower():
            file.write('COORD_FILE_NAME coord-%d.xyz\n'%step)
        elif 'WFN_RESTART_FILE_NAME'.lower() in lines[i].lower():
            file.write('WFN_RESTART_FILE_NAME OT_%d-RESTART.wfn\n'%step)
        elif 'SCF_GUESS'.lower() in lines[i].lower():
            file.write('SCF_GUESS RESTART\n')
        elif 'FILENAME'.lower() in lines[i].lower():
            file.write('FILENAME libra\n')
        else:
            file.write(lines[i])

    file.close()


def run_cp2k_xtb(params):
    """
    This function is for running the CP2K for xTB inputs.

    Args:

        params (dictionary): A dictionary containing the following parameters

          * **cp2k_ot_input_template** (string): The CP2K OT input template file name.
 
          * **cp2k_diag_input_template** (string): The CP2K diagonalization input template file name.

          * **trajectory_xyz_filename** (string): The trajectory xyz file name.

          * **step** (integer): The time step.

          * **cp2k_exe** (string): The full path to CP2K executable.

          * **nprocs** (integer): Number of processors.

    """

    print('**************** Running CP2K ****************')
    ot_input_template = params['cp2k_ot_input_template']
    diag_input_template = params['cp2k_diag_input_template']
    trajectory_xyz_filename = params['trajectory_xyz_filename']
    step = params['step']
    cp2k_exe = params['cp2k_exe']
    mpi_executable = params['mpi_executable']
    nprocs = params['nprocs']

    ##### Run OT
    print('Step', step,'Computing the OT method wfn file...')
    cp2k_xtb_ot_inp(ot_input_template, trajectory_xyz_filename, step)
    t1 = time.time()
    os.system('%s -n %d %s -i xtb_ot_step_%d.inp -o OUT-ot_%d.log'%(mpi_executable, nprocs, cp2k_exe, step, step))
    print('Done with OT wfn. Elapsed time:',time.time()-t1)

    ##### Run diagonalization
    t1 = time.time()
    print('Computing the wfn file using diagonalization method...')
    cp2k_xtb_diag_inp(diag_input_template, step)
    os.system('%s -n %d %s -i xtb_diag_step_%d.inp -o OUT-diag_%d.log'%(mpi_executable, nprocs, cp2k_exe, step, step))
    print('Done with diagonalization. Elapsed time:', time.time()-t1)


def distribute_cp2k_libint_jobs(submit_template: str, run_python_file: str, istep: int, fstep: int, njobs: int, run_slurm: bool, submission_exe='sbatch'):
    """
    This function distributes the jobs to perform CP2K calculations and computing and saving the MO overlaps.

    Args:

        submit_template (string): The name of the submit template to submit a job.

        run_python_file (string): The name of the python file for running the run_cp2k_xtb_step2 function which
                                  contains the parameters.

        istep (integer): The initial step of the MD trajectory.

        fstep (integer): The final step of the MD trajectory.

        njobs (integer): The number of jobs.

        run_slurm (bool): The flag for running the computations either as bash or submitting through sbatch.

        submission_exe (string): The submission executable. For slurm envs, it is 'sbatch' and for pbs envs, it is 'qsub'.

    Return:
        None: but performs the action
    """

    file = open(run_python_file,'r')
    lines = file.readlines()
    file.close()

    nsteps_job = int((fstep-istep)/njobs)
    for njob in range(njobs):
        istep_job = njob*nsteps_job+istep
        fstep_job = (njob+1)*nsteps_job+istep+1
        if njob==(njobs-1):
            fstep_job = fstep

        print('Submitting job',njob+1)
        print('Job',njob,'istep',istep_job,'fstep',fstep_job,'nsteps',fstep_job-istep_job)

        if os.path.exists('job%d'%(njob+1)):
            os.system('rm -rf job%d'%(njob+1))
        os.system('mkdir job%d'%(njob+1))
        os.chdir('job%d'%(njob+1))
        os.system('cp ../%s %s'%(submit_template, submit_template))
        file = open('run.py','w')

        for i in range(len(lines)):
            if 'istep' in lines[i]:
                file.write("params['istep'] = %d\n"%istep_job)
            elif 'fstep' in lines[i]:
                file.write("params['fstep'] = %d\n"%fstep_job)
            else:
                file.write(lines[i])
        file.close()
        if run_slurm:
            os.system('%s %s'%(submission_exe,submit_template))
        else:
            # Just in case you want to use a bash file and not submitting
            os.system('sh %s'%submit_template)
        print('Submitted job', njob)
        os.chdir('../')
        # run_cp2k_xtb_step2(params)



def cp2k_find_excitation_energies(filename):
    """
    This function finds and extracts excitation energies and oscillator strengths 
    from the .log files of the TD-DFPT calculations in CP2k
    
    Args:
        file ( string ): name of the file to read
        
    Returns:
        (list, list):
        
           * the first list contains the energies of excited states
           * the second list contains the oscillator strengths
           
    """
    
    f = open(filename, "r")
    A = f.readlines()
    f.close()
    
    E, F = [], []
    for a in A:
        tmp = a.split()
        
        if len(tmp)== 7:
            if tmp[0] == "TDDFPT|":
                istate = int(float(tmp[1]))
                e = float(tmp[2])
                f = float(tmp[6])
                
                E.append(e)
                F.append(f)
                                
    return E, F

