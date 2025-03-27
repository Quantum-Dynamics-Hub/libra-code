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
import scipy.sparse as sp
import scipy.linalg
import time
import glob 
from libra_py.workflows.nbra import step2_many_body
import scipy.io as io

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn   

from libra_py import data_outs, data_stat, data_conv
from libra_py import units
from libra_py import molden_methods


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
    default_params = { "tolerance":0.05, "isUKS": 0, "lowest_orbital": 1, "highest_orbital": int(1e10)}
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
            if 'Molecular' not in lines[i-1]:
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
                    if int(tmp_splitted_line[0])>=params["lowest_orbital"] and int(tmp_splitted_line[1])<=params["highest_orbital"]:
                        # We need to remove the paranthesis from the 2nd element of the temporary splitted line
                        tmp_spin.append( tmp_splitted_line[1].replace('(','').replace(')','') )
                        tmp_state.append( [ int( tmp_splitted_line[0] ), int( tmp_splitted_line[2] ) ]  )
                        tmp_state_coefficients.append( ci_coefficient  )

            # Here, we have the spin-unpolarize Kohn-Sham basis
            # For this case, spin-components will just return all alpha
            else:
                ci_coefficient = float( tmp_splitted_line[2] )
                if ci_coefficient**2 > tolerance:
                    if int(tmp_splitted_line[0])>=params["lowest_orbital"] and int(tmp_splitted_line[1])<=params["highest_orbital"]:
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
    
        (list, MATRIX(ndof, 1)): labels of all atoms, and their coordinates, ndof = 3 * natoms
	
    """

    f = open(file_name,'r')
    lines = f.readlines()
    f.close()

    # The number of atoms for each time step in the .xyz file of the trajectory.
    number_of_atoms = int(lines[0].split()[0])

    q = MATRIX(3*number_of_atoms, 1)

    # Write the coordinates of the 't' th step in file coord-t.xyz
    f = open('coord-%d'%step+'.xyz','w')

    # This is used to skip the first two lines for each time step.
    n = number_of_atoms+2

    # Write the coordinates of the 'step'th time step into the file
    for i in range( n * step, n * ( step + 1 ) ):
        f.write( lines[i] )
    f.close()

    labels = []
    for i in range(number_of_atoms):
        tmp = lines[ n * step + 2 + i ].split()
        labels.append( tmp[0])
        x = float( tmp[1])  * units.Angst  # convert Angstrom to Bohr
        y = float( tmp[2])  * units.Angst  # convert Angstrom to Bohr
        z = float( tmp[3])  * units.Angst  # convert Angstrom to Bohr
        q.set(3*i + 0, 0,  x)
        q.set(3*i + 1, 0,  y)
        q.set(3*i + 2, 0,  z)

    return labels, q



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


def read_homo_index(filename: str, isUKS: bool = False):
    """
    This function extract the HOMO index from CP2K log file. This index starts from 1.
    Args:

        filename (string): The full path to the log file.

        isUKS (bool): Return 2 indexes for alpha and beta.

    Returns:

        homos (list): The HOMO indexes (the number of occupied orbitals) in the form of [alpha_homo, beta_homo]
        or
        homo (int) The HOMO index (in the case of spin-restricted calculations)
    """
    homos = []
    with open(filename, 'r') as f:
        for line in f:
            if 'occupied' in line.lower() and 'number' in line.lower():
                homos.append(int(line.split()[-1]))
            if len(homos) >= 2:
                break
    if isUKS:
        return homos
    return homos[0]



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
    properties = []
    data = []
    for i in range(len(mat_lines)):
        print(lines[mat_lines[i]])
        properties.append(lines[mat_lines[i]].split()[0])
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
        
    return properties, data


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

    #nsteps_job = int((fstep-istep)/njobs)
    nsteps_job = np.linspace(istep, fstep, njobs+1, endpoint=True, dtype=int).tolist()
    for njob in range(njobs):
        #istep_job = njob*nsteps_job+istep
        #fstep_job = (njob+1)*nsteps_job+istep+1
        istep_job = nsteps_job[njob]
        fstep_job = nsteps_job[njob+1]+1
        #if njob==(njobs-1):
        #    fstep_job = fstep

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



def gaussian_function(a, mu, sigma, num_points, x_min, x_max):
    """
    This function computes the value of a Gaussian function for a set of points.
    Args:
        a (float): The prefactor
        mu (float): The average value, \mu
        sigma (float): The standard deviation, \sigma
        num_points (integer): The number of points in the grid 
        x_min (float): The minimum value of x in the grid
        x_max (float): The maximum value of x in the grid
    Returns:
        x (nparray): The grid
        gaussian_fun (nparray): The computed values of Gaussian function in the grid points, x.
    """
    pre_fact = (a/sigma)/(np.sqrt(2*np.pi))
    x = np.linspace(x_min, x_max, num_points)
    x_input = np.array((-1/2)/(np.square(sigma))*np.square(x-mu))
    gaussian_fun = pre_fact*np.exp(x_input)
    
    return x, gaussian_fun
    
def gaussian_function_vector(a_vec, mu_vec, sigma, num_points, x_min, x_max):
    """
    This function computes the Gaussian function for a set of grid sets using `gaussian_function` and sums them together.
    Args:
        a_vec (list): A list containing the values of prefactors
        mu_vec (list): A list containing the values of averages
        sigma (float): The standard deviation 
        num_points (integer): The number of points in the grid 
        x_min (float): The minimum value of x in the grid
        x_max (float): The maximum value of x in the grid
    Returns:
        x (nparray): The grid
        sum_vec (nparray): The summed values of the Gaussian functions for grids.
    """
    for i in range(len(a_vec)):
        if i==0:
            sum_vec = np.zeros(num_points)
        x, conv_vec = gaussian_function(a_vec[i], mu_vec[i], sigma, num_points, x_min, x_max)
        sum_vec += conv_vec
    return x, sum_vec



def aux_pdos(c1, atoms, pdos_files, margin, homo_occ, orbitals_cols, sigma, npoints, ave_pdos_convolved_all, orbitals, labels):
    """
    c1 - index of the atom kinds
    atoms - the 
    pdos_files (list of strings): names of the files containing pdos data
    shift (double): the amount of "padding" around the minimal and maximal values of the energy scale [units: eV]
    homo_occ (double): occupancy of the HOMO: 2 - for non-polarized, 1 - for polarized
    orbital_cols (list of lists of integers): integers that define which columns correspond to which types of orbitals
    sigma (double): broadening of the pDOS
    npoints (int): how many points to use to compte the convolved pDOS
 
    """

    pdos_ave = np.zeros(np.loadtxt(pdos_files[0]).shape)

    cnt = len(pdos_files)
    for c2, pdos_file in enumerate(pdos_files):
        pdos_mat = np.loadtxt(pdos_file)
        if c2==0:
            pdos_ave = np.zeros(pdos_mat.shape)
        pdos_ave += pdos_mat
    pdos_ave /= cnt

    pdos_ave[:,1] *= units.au2ev
    e_min = np.min(pdos_ave[:,1]) - margin
    e_max = np.max(pdos_ave[:,1]) + margin
    homo_level = np.max(np.where(pdos_ave[:,2]==homo_occ))
    homo_energy = pdos_ave[:,1][homo_level]
    #labels = []

    for c3, orbital_cols in enumerate(orbitals_cols):
        try:
            sum_pdos_ave = np.sum(pdos_ave[:,orbital_cols],axis=1)
            ave_energy_grid, ave_pdos_convolved = gaussian_function_vector(sum_pdos_ave, pdos_ave[:,1], sigma,npoints, e_min, e_max)
            ave_pdos_convolved_all.append(ave_pdos_convolved)
            pdos_label = F"{atoms[1][c1]}, {orbitals[c3]}"
            labels.append(pdos_label)
        except:
            pass

    return ave_energy_grid, homo_energy, ave_pdos_convolved_all, labels


def pdos(params):
    """
    This function computes the weighted pdos for a set of CP2K pdos files.
    Args:
        params (dict):
            * path_to_all_pdos (str): The full path to the pdos files produced by CP2K
            * atoms (list): This list contains the atom types number and names in a separate list.
                            CP2K produces the pdos files in *k{i}*.pdos files where `i` is the atom type number.
                            e.g.: [ [1,2] , ['Ti', 'O'] ]
            * orbitals (list): A list containing the atomic angular momentum.
            * orbital_cols (list): The columns for specific angular momentum as in `orbitals` 
                            e.g. For orbitals of ['s', 'p',        'd',         'f']
                                 [                [3], range(4,7), range(7,12), range(12,19)  ]
            * npoints (integer): The number of grid points for convolution with Gaussian functions.
            * sigma (float): The standard deviation value
            * shift (float): The amount of shifting the grid from the minimum and maximum energy values.
            * is_spin_polarized (int): 0 - non-spin-polarized, 1 - spin-polarized
    Returns:
        ave_energy_grid (nparray): The average energy grid for all pdos files
        homo_energy (float): The value of HOMO energy in eV
        ave_pdos_convolved_all (li st): The list containing all the values 
        labels (lis t): A list containing the labels of the orbital resolved pdos
                       `ave_pdos_convolved_all` is sorted based on the labels in this list
        ave_pdos_convolved_total (nparray): The total density of states
    """
    # Critical parameters
    critical_params = [ "path_to_all_pdos", "atoms"]
    # Default parameters
    default_params = { "orbitals_cols": [[3], range(4,7), range(7,12), range(12,19)], 
                       "orbitals": ['s', 'p', 'd', 'f'], 
                       "sigma": 0.05, "shift": 2.0, "npoints": 2000, "is_spin_polarized": 0
                     } 
    # Check input
    comn.check_input(params, default_params, critical_params) 

    path_to_all_pdos = params['path_to_all_pdos']
    atoms = params['atoms']
    orbitals_cols = params['orbitals_cols']
    orbitals = params['orbitals']
    npoints = params['npoints']
    sigma = params['sigma']
    shift = params['shift']
    is_spin_polarized = params["is_spin_polarized"]


    if is_spin_polarized == 0: # non-polarized case
        homo_occ = 2.0
        labels = []
        ave_pdos_convolved_all = []

        for c1,i in enumerate(atoms[0]):
            # Finding all the pdos files
            pdos_files = glob.glob(path_to_all_pdos+F'/*k{i}*.pdos')

            ave_energy_grid, homo_energy, ave_pdos_convolved_all, labels = aux_pdos(c1, atoms, pdos_files, shift, homo_occ, orbitals_cols, sigma, npoints, ave_pdos_convolved_all, orbitals, labels)

            """
            # Finding all the pdos files
            pdos_files = glob.glob(path_to_all_pdos+F'/*k{i}*.pdos')
            pdos_ave = np.zeros(np.loadtxt(pdos_files[0]).shape)

            cnt = len(pdos_files)
            for c2, pdos_file in enumerate(pdos_files):
                pdos_mat = np.loadtxt(pdos_file)
                if c2==0:
                    pdos_ave = np.zeros(pdos_mat.shape)
                pdos_ave += pdos_mat
            pdos_ave /= cnt
            pdos_ave[:,1] *= units.au2ev
            e_min = np.min(pdos_ave[:,1])-shift
            e_max = np.max(pdos_ave[:,1])+shift
            homo_level = np.max(np.where(pdos_ave[:,2]==2.0))
            homo_energy = pdos_ave[:,1][homo_level]
            for c3, orbital_cols in enumerate(orbitals_cols):
                try:
                    sum_pdos_ave = np.sum(pdos_ave[:,orbital_cols],axis=1)
                    ave_energy_grid, ave_pdos_convolved = gaussian_function_vector(sum_pdos_ave, pdos_ave[:,1], sigma,
                                                                                   npoints, e_min, e_max)
                    ave_pdos_convolved_all.append(ave_pdos_convolved)
                    pdos_label = atoms[1][c1]+F', {orbitals[c3]}'
                    labels.append(pdos_label)
                except:
                    pass
            """

        ave_pdos_convolved_all = np.array(ave_pdos_convolved_all) 
        ave_pdos_convolved_total = np.sum(ave_pdos_convolved_all, axis=0)

        return ave_energy_grid, homo_energy, ave_pdos_convolved_all, labels, ave_pdos_convolved_total

    else:  # spin-polarized case
        homo_occ = 1.0

        labels_alp, labels_bet = [], []
        ave_pdos_convolved_all_alp, ave_pdos_convolved_all_bet = [], []

        for c1,i in enumerate(atoms[0]):
            # Finding all the pdos files
            pdos_files_alp = glob.glob(path_to_all_pdos+F'/*ALPHA*k{i}*.pdos')
            pdos_files_bet = glob.glob(path_to_all_pdos+F'/*BETA*k{i}*.pdos')

            ave_energy_grid_alp, homo_energy_alp, ave_pdos_convolved_all_alp, labels_alp = aux_pdos(c1, atoms, pdos_files_alp, shift, homo_occ, orbitals_cols, sigma, npoints, ave_pdos_convolved_all_alp, orbitals, labels_alp)

            ave_energy_grid_bet, homo_energy_bet, ave_pdos_convolved_all_bet, labels_bet = aux_pdos(c1, atoms, pdos_files_bet, shift, homo_occ, orbitals_cols, sigma, npoints, ave_pdos_convolved_all_bet, orbitals, labels_bet)


        ave_pdos_convolved_all_alp = np.array(ave_pdos_convolved_all_alp)
        ave_pdos_convolved_total_alp = np.sum(ave_pdos_convolved_all_alp, axis=0)

        ave_pdos_convolved_all_bet = np.array(ave_pdos_convolved_all_bet)
        ave_pdos_convolved_total_bet = np.sum(ave_pdos_convolved_all_bet, axis=0)

        return ave_energy_grid_alp, homo_energy_alp, ave_pdos_convolved_all_alp, labels_alp, ave_pdos_convolved_total_alp, \
               ave_energy_grid_bet, homo_energy_bet, ave_pdos_convolved_all_bet, labels_bet, ave_pdos_convolved_total_bet



def exc_analysis(params):
    """
    This function computes the average configuration interaction coefficients for excited states (originally written by Brendan Smith).
    Args:
        params (dict): 
            * path_to_logfiles (string): The full path to all logfiles
            * number_of_states (integer): The number of excited states appeared in the excitation analysis of CP2K outputs
            * tolerance (float): the tolerance value for CI coefficients
            * nsds (integer): The number of SDs to be considered for plotting
            * isUKS (integer): A flag for considering restricted or unrestricted spin calculations:
                               - 0 : No UKS
                               - 1 : With UKS
    Returns:
        coeffs_avg (list): The average CI coefficients
        coeffs_error (list): The error bars for average CI coefficients
    """
    # Critical parameters
    critical_params = [ "path_to_logfiles", "number_of_states"]
    # Default parameters
    default_params = { "isUKS":0, "tolerance": 0.01, "nsds": 2}
    # Check input
    comn.check_input(params, default_params, critical_params)

    path_to_logfiles = params['path_to_logfiles']
    nsds = params['nsds']
    nstates = params['number_of_states']
    
    logfiles = glob.glob(F"{path_to_logfiles}/*.log")
    params.update({'logfile_name': ''})

    ci_coeffs = []
    for logfile in logfiles:
        params.update({"logfile_name": logfile})
        excitation_energies, ci_basis_raw, ci_coefficients_raw_unnorm, spin_components = read_cp2k_tddfpt_log_file( params ) 
        ci_coefficients_raw_norm = step2_many_body.normalize_ci_coefficients(ci_coefficients_raw_unnorm)
        for j in range(len(ci_coefficients_raw_norm)):
            for k in range(len(ci_coefficients_raw_norm[j])):
                ci_coefficients_raw_norm[j][k] = ci_coefficients_raw_norm[j][k]**2
        ci_coeffs.append(ci_coefficients_raw_norm)
    
    
    nsteps = len(ci_coeffs)
    
    coeffs = []
    coeffs_avg   = []
    coeffs_error = []
    
    
    for state in range(nstates):
        coeffs.append( [] )
        coeffs_avg.append( [] )
        coeffs_error.append( [] )
    
        for sd in range( nsds ):
    
            coeffs[state].append( [] )
            coeffs_avg[state].append( [] )
            coeffs_error[state].append( [] )
    
            for step in range( nsteps ):
                if len( ci_coeffs[step][state] ) < nsds and sd > len( ci_coeffs[step][state] )-1:
                    coeffs[state][sd].append( 0.0 )
                else:
                    coeffs[state][sd].append( ci_coeffs[step][state][sd] )
         
            mb_coeff_avg, mb_coeff_std = data_stat.scalar_stat( coeffs[state][sd] )
            coeffs_avg[state][sd].append( mb_coeff_avg )
            coeffs_error[state][sd].append( 1.96 * mb_coeff_std / np.sqrt(nsteps) )
      

    return coeffs_avg, coeffs_error


def extract_energies_sparse(params):
    """
    This function extract the energies vs time data from scipy.sparse files
    Args:
        params (dict):
            * path_to_energy_files (string): The full path to sparse energy files 
            * dt (float): The time step
            * prefix (string): The prefix of the files
            * suffix (string): The suffix of the files
            * istep (integer): The initial step
            * fstep (integer): The final step
    Returns:
        md_time (nparray): Contains the md time and can be used for plotting
        energies (nparray): The energies of each state
    """
    critical_params = [ "path_to_energy_files"]
    # Default parameters
    default_params = { "dt":1.0, "prefix": "Hvib_sd_", "suffix": "_re", "istep": 0, "fstep": 5}
    # Check input
    comn.check_input(params, default_params, critical_params)

    path_to_energy_files = params['path_to_energy_files']
    dt = params['dt']
    prefix = params['prefix']
    suffix = params['suffix']
    istep = params['istep']
    fstep = params['fstep']
    
    energies = []
    for step in range(istep,fstep):
        file = F'{path_to_energy_files}/{prefix}{step}{suffix}.npz'
        energies.append(np.diag(sp.load_npz(file).todense().real))
    energies = np.array(energies)
    md_time = np.arange(0,energies.shape[0]*dt,dt)
    
    return md_time, energies


def read_wfn_file(filename):
    """
    This function reads the data stored in wfn files produced by CP2K.
    These data are stored in this format (from the cp2k/src/qs_mo_io.F file):
    
    From the subroutine 'write_mo_set_low':
    ### ================================================== 
    WRITE (ires) natom, nspin, nao, nset_max, nshell_max
    WRITE (ires) nset_info
    WRITE (ires) nshell_info
    WRITE (ires) nso_info
    
    DO ispin = 1, nspin
         nmo = mo_array(ispin)%nmo
         IF ((ires > 0) .AND. (nmo > 0)) THEN
            WRITE (ires) nmo, &
               mo_array(ispin)%homo, &
               mo_array(ispin)%lfomo, &
               mo_array(ispin)%nelectron
            WRITE (ires) mo_array(ispin)%eigenvalues(1:nmo), &
               mo_array(ispin)%occupation_numbers(1:nmo)
         END IF
         IF (PRESENT(rt_mos)) THEN
            DO imat = 2*ispin - 1, 2*ispin
               CALL cp_fm_write_unformatted(rt_mos(imat), ires)
            END DO
         ELSE
            CALL cp_fm_write_unformatted(mo_array(ispin)%mo_coeff, ires)
         END IF
      END DO
    ### ================================================== 
    SUMMARY:
    natom, nspin, nao, nset_max, nshell_max, nset_info, nshell_info, nso_info
    for ispin in range(nspins):
        nmo, homo, lfomo, nelectron, (eigenvalues(1:nmo), occupation_numbers(1:nmo)), mo_coeff

    Args:
        filename (string): The name of the wfn file

    Returns:
        basis_data (list): A list that contains the basis set data
        spin_data (list): A list containing the alpha and beta spin data
        eigen_vals_and_occ_nums (list): A list that contains the energies and occupation numbers.
    """
    wfn_file = io.FortranFile(filename,'r')
    basis_data = []
    natom, nspin, nao, nset_max, nshell_max = wfn_file.read_ints()
    nset_info = wfn_file.read_ints()
    nshell_info = wfn_file.read_ints()
    nso_info = wfn_file.read_ints()
    basis_data.append(natom)
    basis_data.append(nspin)
    basis_data.append(nao)
    basis_data.append(nset_max)
    basis_data.append(nshell_max)
    basis_data.append(nset_info)
    basis_data.append(nshell_info)
    basis_data.append(nso_info)
    spin_data = []
    eigen_vals_and_occ_nums = []
    mo_coeffs = []
    for spin in range(nspin):
        spin_tmp = []
        nmo, homo, lfomo, nelectron = wfn_file.read_ints()
        spin_tmp.append(nmo)
        spin_tmp.append(homo)
        spin_tmp.append(lfomo)
        spin_tmp.append(nelectron)
        spin_data.append(spin_tmp)
        eigen_val_and_occ_num = wfn_file.read_reals(dtype=np.float32)
        eigen_vals_and_occ_nums.append(eigen_val_and_occ_num)
        mo_coeff = []
        for mo in range(nmo):
            coeffs = wfn_file.read_reals()
            mo_coeff.append(coeffs)
        mo_coeffs.append(mo_coeff)
    wfn_file.close()
    return basis_data, spin_data, eigen_vals_and_occ_nums, mo_coeffs


def write_wfn_file(filename, basis_data, spin_data, eigen_vals_and_occ_nums, mo_coeffs):
    """
    This function writes a set of data including basis sets, spin data, eigenvalues and occupation numbers, and 
    MO coefficients into a binary file (wfn) readable by CP2K.
    Args:
        filename (string): The name of the output wfn file.
        basis_data (list): A list that contains the basis set data
        spin_data (list): A list containing the alpha and beta spin data
        eigen_vals_and_occ_nums (list): A list that contains the energies and occupation numbers.
        mo_coeffs (numpy array): The numpy array that contains the MO coefficients.
    Returns:
        None
    """
    wfn_file = io.FortranFile(filename,'w')
    natom = basis_data[0]
    nspin = basis_data[1]
    nao = basis_data[2]
    nset_max = basis_data[3]
    nshell_max = basis_data[4]
    nset_info = basis_data[5]
    nshell_info = basis_data[6]
    nso_info = basis_data[7]
    wfn_file.write_record(natom, nspin, nao, nset_max, nshell_max)
    wfn_file.write_record(nset_info)
    wfn_file.write_record(nshell_info)
    wfn_file.write_record(nso_info)
    for spin in range(nspin):
        wfn_file.write_record(spin_data[spin][0], spin_data[spin][1], spin_data[spin][2], spin_data[spin][3])
        nmo = spin_data[spin][0]
        wfn_file.write_record(eigen_vals_and_occ_nums[spin])#, dtype=np.float32)
        #zero_eig_occ = np.zeros(eigen_vals_and_occ_nums[spin].shape)
        #wfn_file.write_record(zero_eig_occ)#, dtype=np.float32)
        for mo in range(nmo):
            wfn_file.write_record(mo_coeffs[spin][mo])
            #zero_coeffs = np.random.random(mo_coeffs[spin][mo].shape)
            #wfn_file.write_record(zero_coeffs)
    wfn_file.close()


def compute_energies_coeffs(ks_mat, overlap):
    """
    This function solves the general eigenvalue problem described above using a Cholesky decomposition
    of the overlap matrix. The eigenvalues are sorted.
    More information: https://doi.org/10.1016/j.cpc.2004.12.014
    Args:
        ks_mat (numpy array): The Kohn-Sham matrix
        overlap (numpy array): The atomic orbital overlap matrix
    Returns:
        eigenvalues (numpy array): The energies (eigenvalues)
        eigenvectors (numpy array): The MO coefficients
    """
    # Cholesky decomposition of the overlap matrix
    U = np.linalg.cholesky( overlap ).T
    # One ca also use the following as well but it is computationally more demanding
    # U = scipy.linalg.fractional_matrix_power(S, 0.5)
    U_inv = np.linalg.inv( U )
    UT_inv = np.linalg.inv( U.T )
    K_prime = np.linalg.multi_dot( [UT_inv, ks_mat, U_inv] )
    eigenvalues, eigenvectors = np.linalg.eig( K_prime )
    # Transform back the coefficients 
    eigenvectors = np.dot(U_inv, eigenvectors)
    sorted_indices = np.argsort(eigenvalues)
    eigenvectors = eigenvectors[:,sorted_indices].T
    
    
    return eigenvalues[sorted_indices], eigenvectors


def compute_energies_coeffs_scipy(ks_mat, overlap):
    """
    This function solves the general eigenvalue problem described above using a Cholesky decomposition
    of the overlap matrix. The eigenvalues are sorted.
    More information: https://doi.org/10.1016/j.cpc.2004.12.014
    Args:
        ks_mat (numpy array): The Kohn-Sham matrix
        overlap (numpy array): The atomic orbital overlap matrix
    Returns:
        eigenvalues (numpy array): The energies (eigenvalues)
        eigenvectors (numpy array): The MO coefficients
    """
    # Cholesky decomposition of the overlap matrix
    U = scipy.linalg.cholesky( overlap ).T
    # One ca also use the following as well but it is computationally more demanding
    # U = scipy.linalg.fractional_matrix_power(S, 0.5)
    U_inv = scipy.linalg.inv( U )
    UT_inv = scipy.linalg.inv( U.T )
    #K_prime = scipy.linalg.multi_dot( [UT_inv, ks_mat, U_inv] )
    K_prime = UT_inv @ ks_mat @ U_inv
    eigenvalues, eigenvectors = scipy.linalg.eig( K_prime )
    # Transform back the coefficients 
    eigenvectors = U_inv @ eigenvectors
    sorted_indices = np.argsort(eigenvalues) 
    eigenvectors = eigenvectors[:,sorted_indices].T
    
    
    return eigenvalues[sorted_indices], eigenvectors

def compute_density_matrix(eigenvectors, homo_index):
    """
    This function computes the density matrix: P=2\times c_{occ}\times c_{occ}^T
    Args:
        eigenvectors (numpy array): The set of eigenvectors
        homo_index (integer): The HOMO index to select only the occupied orbitals
    Returns:
        density_mat (numpy array): The density matrix.
    """
    occupied_eigenvectors = eigenvectors[0:homo_index]
    density_mat = 2*np.dot(occupied_eigenvectors.T, occupied_eigenvectors)
    return density_mat

def compute_convergence(ks_mat, overlap, density_mat):
    """
    This function computes the convergence error based on the commutator relation: e=KPS-SPK 
    More information: https://doi.org/10.1016/j.cpc.2004.12.014
    At this point, this function gives not exact results as in CP2K and the results seems to be different.
    TODO: Need to check this in more details by checking the CP2K source code.
    """
    return np.linalg.multi_dot([ks_mat,density_mat,overlap])-np.linalg.multi_dot([overlap,density_mat,ks_mat])

def scf_timing(log_filename):
    """
    This function returns the timing data related to CP2K log file.
    Args:
        log_filename (string): The name of the logfile
    Returns:
        ncycle (int): Number of SCF cycles in the log file
        convergence (numpy array): The convergence value at each SCF step 
        timing (float): The total time of the SCF cycle
    """
    f = open(log_filename,'r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
        if '*** SCF run ' in lines[i] or 'Leaving inner SCF loop' in lines[i]:
            end = i
        elif 'Step' in lines[i] and 'Convergence' in lines[i]:
            start = i+2
    timing = 0
    ncycle = 0
    convergence = []
    for i in range(start, end):
        tmp = lines[i].split()
        if len(tmp)>0:
            ncycle += 1
            timing += float(tmp[3])
            convergence.append(float(tmp[4]))
    return ncycle, np.array(convergence), timing


def resort_ao_matrices(l_vals):
    """
    This function returns the resotring indices for resoting the atomic orbital
    matrices from Psi4 to CP2K order according to this order:
    
    AO matrices order (example for Cd atom):
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
        reordered_indices = index_reorder_v2(l_val)
        for j in range(len(reordered_indices)):
            # now append it by plus the counter since
            # we aim to do it for all the eigenvectors
            # and l values
            new_indices.append(c+reordered_indices[j])
        # increase the counter with t
        c += len(reordered_indices)
    # Return the new indices
    return new_indices
    
    
def index_reorder_v2(l_val):
    """
    This function returns the new ordering of the angular momentum value based on the 
    order used in Psi4. Detailed explanation was given for the resort_ao_matrices function.
    
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
        new_order = [3,1,2]
    # for d orbital
    elif l_val == 2:
        new_order = [5,3,1,2,4]
    # for f orbital
    elif l_val == 3:
        new_order = [7,5,3,1,2,4,6]
    # for g orbital
    elif l_val == 4:
        new_order = [9,7,5,3,1,2,4,6,8]

    # The indeices
    return np.array(new_order)-1



def atom_components_cp2k(filename):
    """
    This function finds the components of the angular momentum of each atom
    in the system and returns a list that can be used to find specific blocks
    in the atomic orbital matrices such as Hamiltonians or overlaps.
    Args:
        filename (string): The name of the AO matrices file
    Returns:
        indices (list): The list of indices related to the angular momentum
                        components of each atom in the system.
    """
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    atoms = []
    for i in range(len(lines)):
        tmp = lines[i].split()
        if len(tmp)>4:
            try:
                atom = int(tmp[1])
                if len(atoms)==0:
                    atoms.append(atom)
                else:
                    if atom>=atoms[-1]:
                        atoms.append(atom)
                    else:
                        break
            except:
                pass
    np.unique(atoms).shape
    indices = []
    for i in np.unique(atoms):
        index = list(np.where(atoms==i)[0])
        indices.append(index)
    return indices

def compute_mo_overlap(params):
    """
    This function computes the molecular orbital overlap calculations between
    two geometries: $\langle\Phi_{i,R}\|\Phi_{j,R'}\rangle$. It takes a dictionary as parameters.
    Args:
        params (dictionary):
            molden_file_1 (string): The path to the first geometry molden file
            molden_file_2 (string): The path to the second geometry molden file
            is_spherical (bool): The spherical coordinate boolean flga
            nprocs (integer): The number of pocessors
    Returns:
        mo_overlap_matrix (numpy array): The overlap between molecular orbitals
    """
    # Create Libint integral shells for each molden file
    shell_1, l_vals_1 = molden_methods.molden_file_to_libint_shell(params['molden_file_1'], params['is_spherical'])
    shell_2, l_vals_2 = molden_methods.molden_file_to_libint_shell(params['molden_file_2'], params['is_spherical'])
    # Reindexing the eigenvectors outputted by CP2K
    new_indices = resort_molog_eigenvectors(l_vals_1)
    eig_vect_1, energies_1 = molden_methods.eigenvectors_molden(params['molden_file_1'], nbasis(shell_1), l_vals_1)
    eig_vect_2, energies_2 = molden_methods.eigenvectors_molden(params['molden_file_2'], nbasis(shell_2), l_vals_2)
    
    # The new eigenvectors
    resortted_eig_vect_1 = []
    resortted_eig_vect_2 = []
    for j in range(len(eig_vect_1)):
        # the new and sorted eigenvector
        eigenvector1 = eig_vect_1[j]
        eigenvector1 = eigenvector1[new_indices]
        eigenvector2 = eig_vect_2[j]
        eigenvector2 = eigenvector2[new_indices]
        # append it to the eigenvectors list
        resortted_eig_vect_1.append(eigenvector1)
        resortted_eig_vect_2.append(eigenvector2)
    resortted_eig_vect_1 = np.array(resortted_eig_vect_1)
    resortted_eig_vect_2 = np.array(resortted_eig_vect_2)
    ao_matrix = compute_overlaps(shell_1, shell_2, params['nprocs'])
    ao_matrix = data_conv.MATRIX2nparray(ao_matrix)
    mo_overlap_matrix = np.linalg.multi_dot([resortted_eig_vect_1, ao_matrix, 
                                             resortted_eig_vect_2.T]) 
    return mo_overlap_matrix

