#*********************************************************************************
#* Copyright (C) 2020-2023 Mohammad Shakiba, Brendan Smith, Alexey V. Akimov
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: Gaussian_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing the Gaussian inputs and outputs.
.. moduleauthors:: 
       Mohammad Shakiba, Brendan Smith, Alexey V. Akimov 
  
"""



import os
import sys
import numpy as np
import libra_py.packages.cp2k.methods as CP2K_methods
import util.libutil as comn


def read_gaussian_tddft_log_file(params):
    """
    This function read Gaussian output file and extracts the excitation analysis results.

    Args:

        params (dictionary): The dictionary containing the input parameters.

            logfile_name (string): The log file name.

            tolerance (float): The tolerance factor for choosing the excited states of the excitation analysis.

            number_of_states (integer): The number of excited states to be considered in the excitation analysis.

            isUKS (integer): This parameter is the flag for restricted or unrestricted spin calculations. If it is
                             set to 1 the unrestricted spin calculations will be considered.
    Returns:

        excitation_energies (1D numpy arra): The excitation energies.

        ci_basis (list): The list containnig the CI basis of the excitation analysis.

        ci_coefficients (list): The list containing the ci_coefficients for each excitation.

        spin_components (list): The list containing the spin components for each excitation.

    """

    # Critical parmeters
    critical_params = [ "logfile_name" ]
    # Default parameters
    default_params = { "tolerance": 0.05, "number_of_states": 5, "isUKS": 0 }
    # Check input
    comn.check_input(params, default_params, critical_params)

    # Gaussian log file name
    gaussian_log_file_name = params["logfile_name"]
    # tolerance factor
    tolerance = params["tolerance"]
    # The unrestricted spin calculations flag
    isUKS = params["isUKS"]
    # The number of states
    number_of_states = params["number_of_states"]

    # Open the log file
    file = open(gaussian_log_file_name,'r')
    lines = file.readlines()
    file.close()

    # Initialize the excitation energies list
    excitation_energies = []
    # The lines with excited state in them
    excited_states_lines = []
    for i in range(len(lines)):
        # If 'Excited State' was found in the log file append it in excited_states_lines
        if 'Excited State' in lines[i]:
            excitation_energies.append(float(lines[i].split()[4]))
            excited_states_lines.append(i)
    print("excited_states_lines = ", excited_states_lines)

    # Initialize the ci_basis, ci_coefficients, and spin_components
    ci_basis = []
    ci_coefficients = []
    spin_components = []
    # Find the final line of the excitation analysis (the blank line) and append
    # it to excited_states_lines for further use in the 'for' loops
    for i in range(excited_states_lines[-1],len(lines)):
        # Find the blank line
        if len(lines[i].split())==0:
            excited_states_lines.append(i)
            # If founf the blank line exit the for loop
            break
    
    for i in range(number_of_states):
        # Initialize tmp variables for storing the excitation analyses
        tmp_ci_state              = []
        tmp_ci_state_coefficients = []
        tmp_spin = []
        
        for j in range( excited_states_lines[i], excited_states_lines[i+1] ):
            if '->' in lines[j]: # and len(lines[j].split())==4:
                if isUKS==1:
                    # For alpha spin
                    if 'A' in lines[j]:
                        line_alpha = lines[j].replace('A','').replace('->','').split()
                        # Use the tolerance factor for chosing the states
                        if float(line_alpha[2])**2 > tolerance:
                            tmp_ci_state.append( [int(line_alpha[0]), int(line_alpha[1])] )
                            tmp_ci_state_coefficients.append( float(line_alpha[2]) )
                            tmp_spin.append('alp')

                    # For beta spin
                    elif 'B' in lines[j]:
                        line_beta = lines[j].replace('B','').replace('->','').split()
                        # Use the tolerance factor for chosing the states
                        if float(line_beta[2])**2 > tolerance:
                            tmp_ci_state.append( [int(line_beta[0]), int(line_beta[1])] )
                            tmp_ci_state_coefficients.append( float(line_beta[2]) )
                            tmp_spin.append('bet')
                else:
                    # Just for alpha spin and the same as above
                    print("\nThis is the spin-restricted case, we are reading the lines now")
                    tmp_line_alpha = lines[j].replace("->","")
                    line_alpha = tmp_line_alpha.split()

                    print("\nlines[j] = ", lines[j])  
                    print("tmp_line_alpha = ", tmp_line_alpha)
                    print("line_alpha = ", line_alpha)

                    if float(line_alpha[2])**2 > tolerance:
                        tmp_ci_state.append( [int(line_alpha[0]), int(line_alpha[1])] )
                        tmp_ci_state_coefficients.append( float(line_alpha[2]) )
                        tmp_spin.append('alp')            

        # Append the tmp variables into the main variables
        ci_basis.append(tmp_ci_state)
        ci_coefficients.append(tmp_ci_state_coefficients)
        spin_components.append(tmp_spin)
    
    return excitation_energies[0:number_of_states], ci_basis, ci_coefficients, spin_components




def read_energies_from_gaussian_log_file( params ):
    """
    This function read the energies from Gaussian log file and returns the Kohn-Sham and total_energies.

    Args:

        params (dictionary): The dictionary containing the input parameters.

            logfile_name (string): The log file name.

            min_band (integer): The minimum state number.

            max_band (integer): The maximum state number.

            spin (integer): The spin component. 1 is for alpha and 2 is for beta spin.

    Returns:

        ks_energies (1D numpy array): The Kohn-Sham energies.

        total_energy (float): The total energy of the system.

    """
    # Critical params
    critical_params = [ "logfile_name", "min_band", "max_band" ]
    # Default params
    default_params = { "spin": 1 }
    # Check inputs
    comn.check_input(params, default_params, critical_params)

    # Extracting the input parameters
    gaussian_log_file_name = params["logfile_name"]
    min_band = params["min_band"]
    max_band = params["max_band"]
    spin     = params["spin"]

    file = open(gaussian_log_file_name,'r')
    lines = file.readlines()
    file.close()

    # Initialize the occupied and unoccupied energies
    occupied_energies   = []
    unoccupied_energies = []
    if spin==1:
        # alpha spin
        spin_letter = 'alpha'
    elif spin==2:
        # beta spin
        spin_letter = 'beta'
    
    for i in range(len(lines)):
        # Find the eigenvalues of the occupied molecular orbitals
        if 'eigenvalues' in lines[i].lower() and 'occ' in lines[i].lower() and spin_letter in lines[i].lower():
            for j in range(len(lines[i].split())):
                try:
                    occupied_energies.append(float(lines[i].split()[j]))
                except:
                    pass
        # Find the eigenvalues of the unoccupied molecular orbitals
        if 'eigenvalues' in lines[i].lower() and 'virt' in lines[i].lower() and spin_letter in lines[i].lower():
            for j in range(len(lines[i].split())):
                try:
                    unoccupied_energies.append(float(lines[i].split()[j]))
                except:
                    pass
    # Turn them into numpy arrays
    occupied_energies = np.array(occupied_energies)
    unoccupied_energies = np.array(unoccupied_energies)
    # Concatenate the occupied and unoccpied energies so we can choose from min_band to max_band
    ks_energies = np.concatenate((occupied_energies,unoccupied_energies))

    # Now the total energy
    total_energy = 0
    for i in range(len(lines)):
        # Find the 'SCF Done' in the output, total energy is there
        if 'SCF Done:'.lower() in lines[i].lower():
            for j in range(len(lines[i].split())):
                try:
                    total_energy = float( lines[i].split()[j] )
                    break
                except:
                    pass
            break
    
    return ks_energies[min_band-1:max_band], total_energy




def gaussian_distribute( istep, fstep, nsteps_this_job, trajectory_xyz_file, gaussian_input, curr_job_number ):
    """
    Distributes Gaussian jobs for trivial parallelization 

    Make sure that your Gaussian input file has absolute paths to the following input parameters:
        This parameters should be set in the input template
        chk file ---> e.g.  %chk=/home/username/gaussian_calculations/check_file.chk
        rwf file ---> e.g.  %rwf=/home/username/gaussian_calculations/check_file.rwf

    Args:
        
        istep (integer): The initial time step in the trajectory xyz file.

        fstep (integer): The final time step in the trajectory xyz file.

        nsteps_this_job (integer): The number of steps for this job.

        trajectory_xyz_file (string): The full path to trajectory xyz file.

        gaussian_input (string): The sample Gaussian input template.

        curr_job_number (integer): The current job number.

    Returns:

        None

    """

    # Now we need to distribute the jobs into job batches
    # First, make a working directory where the calculations will take place
    os.chdir("wd")

    nsteps = fstep - istep + 1
    njobs  = int( nsteps / nsteps_this_job ) 

    # Initialize the curr_step to istep
    curr_step = istep

    # Make the job directory folders
    os.system("mkdir job"+str(curr_job_number)+"")
    os.chdir("job"+str(curr_job_number)+"")
    # Copy the trajectory file and input template there
    os.system("cp ../../"+trajectory_xyz_file+" .")
    os.system("cp ../../"+gaussian_input+" .")

    # Now, we need to edit the submit file
    # Now, in jobs folder njob, we should do only a certain number of steps
    for step in range( nsteps_this_job ):

        # Extract the coordinates and write them to a xyz file
        CP2K_methods.read_trajectory_xyz_file( trajectory_xyz_file, curr_step )

        # Now, we need to edit the gaussian_input file by adding the 
        # coordinates to the input file
        tmp = open(gaussian_input)
        A   = tmp.readlines()
        sz  = len(A)
        tmp.close()
        tmp2 = open("step_%d"%curr_step+".gjf","w"); tmp2.close()
        for i in range(sz):

            b = A[i].strip().split()

            if not b:
                continue

            tmp2 = open("step_%d"%curr_step+".gjf","a")
            tmp2.write(A[i])
            tmp2.close()
        # add to the curr_step
        curr_step += 1

    # Go back to the main directory
    os.chdir("../../") 



def read_trajectory_xyz_file_gaussian(trajectory_xyz_file_name: str, time_step: int):
    """
    This function reads the trajectory of a molecular dynamics .xyz file and
    extract the 'time_step' th step then returns it in form of a numpy array.

    Args:

        trajectory_xyz_file_name (string): The trajectory .xyz file name.

        time_step (integer): The desired time to extract its .xyz 
                     coordinates which starts from zero.

    Returns:

        coordinates (2D numpy array): The numpy array containing the coordinates of the time_step.

    """

    f = open(trajectory_xyz_file_name,'r')
    lines = f.readlines()
    f.close()

    # The number of atoms for each time step in the .xyz file of the trajectory.
    number_of_atoms = int(lines[0].split()[0])

    # Add two to the number of atoms, this is the number of 
    # lines for each geometry in the trajectory xyz file
    n = number_of_atoms+2
    # Initialize the coordiantes
    coordinates = []
    # Find the time_step coordinates
    for i in range( n*time_step, n*( time_step + 1 ) ):
        coordinates.append(lines[i].split())
    
    coordinates = np.array(coordinates)
    
    return coordinates


def gaussian_input( project_name, time_step, sample_input, trajectory_xyz_file_name ):
    """
    This function creates a Gaussian input file for the time_step geometry using an input template.

    Args:

        project_name (string): The project name.

        time_step (integer): The time step.

        sample_input (string): The Gaussian sample input.

        trajectory_xyz_file_name (string): The trajectory xyz file name.

    Returns:

        None

    """ 
    # Read the sample input lines
    f = open(sample_input,'r')
    lines = f.readlines()
    f.close()
    
    for i in range(0,len(lines)):
        # Find the chk line for placing the check file
        if '%chk=' in lines[i].lower():
            chk_line = i
        # Find the rwf line for placing the rwf file
        if '%rwf=' in lines[i].lower():
            rwf_line = i
        # Find the line with '#', this line shows the calculation types for Gaussian
        if '#' in lines[i]:
            calc_type_line = i

    # Extract the coordinates from the trajectory xyz file
    coordinates = read_trajectory_xyz_file_gaussian(trajectory_xyz_file_name, time_step)
    # Open a new input file
    f = open('%s-%d.gjf'%(project_name,time_step),'w')
    
    for i in range(0,calc_type_line+5):
        if i==chk_line:
            d = lines[chk_line].split('/')[-1]
            
            f.write('%%chk=%s'%os.getcwd()+'/'+'%s-%d.chk'%(project_name,time_step))
            f.write('\n')
        elif i==rwf_line:
            d = lines[rwf_line].split('/')[-1]
            f.write('%%rwf=%s'%os.getcwd()+'/'+'%s-%d.rwf'%(project_name,time_step))
            f.write('\n')
        else:
            f.write(lines[i])


    for i in range( 2, len(coordinates) ):
        for j in range( 0, len(coordinates[i]) ):
            f.write( coordinates[i][j] )
            f.write('  ')
        f.write('\n')    
    # We need some blank lines at the end of the input file, here we add three
    f.write('\n')
    f.write('\n')
    f.write('\n')
    # Close the file
    f.close()




def cube_generator_gaussian( project_name, time_step, min_band, max_band, nprocs, sample_cube_file, isUKS ):
    """
    This function generates the cube files by first forming the 'fchk' file from 'chk' file.
    Then it will generate the cube files from min_band to max_band the same as CP2K naming.
    This will helps us to use the read_cube and integrate_cube functions easier.

    Args:

        project_name (string): The project_name.

        time_step (integer): The time step.

        min_band (integer): The minimum state number.

        max_band (integer): The maximum state number.

        nprocs (integer): The number of processors used to generate the cubes using 'cubegen'.

        sample_cube_file (str): The path to a sample cube file. This file is used to generate the cube
                                files according the mesh of this file. Therefore the integration will be plausible.

        isUKS (integer): The unrestricted spin calculation flag. 1 is for spin unrestricted calculations.
                         Other numbers are for spin restricted calculations

    Returns:

        None

    """

    # Form the fchk file from the chk file
    os.system('formchk %s-%d.chk %s-%d.fchk'%(project_name, time_step, project_name, time_step))
    # Print the commands...
    print('formchk %s-%d.chk %s-%d.fchk'%(project_name, time_step, project_name, time_step))
    # Generate the sample cube file, if it exists it will not create it again
    if not os.path.isfile(sample_cube_file):
        # Generate the sample cube file with the maximum fineness ---> 100 as in Gaussian website
        os.system('cubegen %d MO=Homo %s-%d.fchk %s 100 h'%(nprocs, project_name, time_step, sample_cube_file))
    # For spin unrestricted
    if isUKS == 1:
        # Generate the names and cube files for each state. Here we use the cube names by CP2K. This will increase the speed of our work.
        for state in range(min_band,max_band+1):
            # State name in CP2K format
            state_name = CP2K_methods.state_num_cp2k(state)
            # Cube file name 
            cube_name = '%s-%d-WFN_%s_1-1_0.cube'%(project_name, time_step, state_name)
            print('Generating cube for state %d'%state)
            # Generate cube files for alpha spin using the 'cubegen'
            os.system('cubegen %d AMO=%d %s-%d.fchk %s -1 h %s'%(nprocs, state, project_name, time_step, cube_name, sample_cube_file))
            print('cubegen %d AMO=%d %s-%d.fchk %s -1 h %s'%(nprocs, state, project_name, time_step, cube_name, sample_cube_file))
            cube_name = '%s-%d-WFN_%s_2-1_0.cube'%(project_name, time_step, state_name)
            print('Generating cube for state %d'%state)
            # Generate cube files for beta spin using the 'cubegen'
            os.system('cubegen %d BMO=%d %s-%d.fchk %s -1 h %s'%(nprocs, state, project_name, time_step, cube_name, sample_cube_file))
            print('cubegen %d BMO=%d %s-%d.fchk %s -1 h %s'%(nprocs, state, project_name, time_step, cube_name, sample_cube_file))
    # Spin restricted case
    else:
        for state in range(min_band,max_band+1):
            # Use cp2k names because the rest of the code expects this format
            state_name = CP2K_methods.state_num_cp2k(state)
            cube_name = '%s-%d-WFN_%s_1-1_0.cube'%(project_name, time_step, state_name)
            print('Generating cube for state %d'%state)
            # Generate the cubes for alpha spin only
            os.system('cubegen %d MO=%d %s-%d.fchk %s -1 h %s'%(nprocs, state, project_name, time_step, cube_name, sample_cube_file))
            print('cubegen %d MO=%d %s-%d.fchk %s -1 h %s'%(nprocs, state, project_name, time_step, cube_name, sample_cube_file))
