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
.. module:: cube_file_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing with cube files
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




def read_cube(filename: str):
    """
    This function reads the wavefunction from a cube file and stores it in
    a 1D numpy array

    Args:

        filename (string): the name of the .cube file to read

    Returns:

        data (numpy.array) : the 1D array of the wavefunctions for all the points on the grid
	
    """
    
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    natoms = int(lines[2].split()[0])
    # We skip a few lines in the cube files, that go as follows:
    # 2 lines - comments
    # 1 line  - the number of atoms, etc.
    # 3 lines - the grid spacing and number of grid points in each dimensions

    nstart = natoms+2+1+3        # the index of the first line containing wfc data
    nlines = len(lines)          # the total number of lines

    x = []
    for i in range(nstart,nlines):
        tmp = lines[i].split()
        ncols = len(tmp)

        for j in range(ncols):
            x.append(float(tmp[j]))
    
    data = np.array(x)
    #data = np.loadtxt(filename,skiprows=n)
    
    return data




def grid_volume(filename: str):
    """
	This function reads the wavefunction from a cube file and calculate 
	the grid volum using the X-, Y- and Z-axis of the volumetric region
	which are placed in the 4th, 5th and 6th line of the cube file structure
	
        Args:
        
	    filename (string): The name of the .cube file to read.
		
        Returns:
	
	    dv (float): The grid volume in Bohr^3.
        	
    """
	
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    
    # We use the 3rd, 4th and 5th row in the lines to obtain
    # the axes of the parallelpiped into a numpy array (Voxel).
    axis_1 = [float(lines[3].split()[1]),float(lines[3].split()[2]),float(lines[3].split()[3])]
    axis_2 = [float(lines[4].split()[1]),float(lines[4].split()[2]),float(lines[4].split()[3])]
    axis_3 = [float(lines[5].split()[1]),float(lines[5].split()[2]),float(lines[5].split()[3])]
    vol_element = np.array([axis_1, axis_2, axis_3])
    # Then we calculate the determinant of Voxel to obtain the volume.
    dv = np.absolute(np.linalg.det(vol_element))
	
    
    return dv 




def read_volumetric_data(filename: str):
    """
    This function reads the volumetric data in a format used for plotting
    of the data. The difference between 'read_cube' function is that it will
    show the data in a 3D array which has the shape as the number of grid points
    for each of the X-, Y-, and Z-axis in the 4th to 6th line of the cube file.
	
    This function will return the grid points in each axis, the coordinates of the
    structure and the spacing vector which is used to plot the isosurfaces of the 
    molecular orbitals.
	
    Args:
	
        filename (string): The name of the .cube file.
		
    Returns:
	
        coordinates (2D numpy array): The coordinates of the molecule structure in
                                      the same format shown in the .cube file.
		
        x_grid, y_grid, z_grid (3D numpy array): Containing the grid points for each of the X-,
                                     Y-, and Z- axis.
		
        wave_fun (3D numpy array): The volumetric data in a 3D numpy array format.
		
        
	spacing_vector (numpy 1D array): The spacing vector used for plotting the isosurfaces.

    """

    f = open(filename,'r')
    lines = f.readlines()
    f.close()
	
    # The number of atoms in the 3rd line
    natoms = int(lines[2].split()[0])


    # The number of voxels defined for each axis obtained from
    # the first elements of the 4th, 5th, and 6th line of the cube files
    nx = int(lines[3].split()[0])
    ny = int(lines[4].split()[0])
    nz = int(lines[5].split()[0])


    # The three vectors below are the same vectors which are present in lines 4th to 6th
    # which are used to create the 3D numpy array of the grid points (x, y, z) used to plot
    # the isosurfaces.
    # Here we use the same unit as is used in the .cube file structure which is Bohr
    # with no need for unit conversion. This makes the plotting easier.
    axis_1 = np.array([float(lines[3].split()[1]),float(lines[3].split()[2]),float(lines[3].split()[3])])
    axis_2 = np.array([float(lines[4].split()[1]),float(lines[4].split()[2]),float(lines[4].split()[3])])
    axis_3 = np.array([float(lines[5].split()[1]),float(lines[5].split()[2]),float(lines[5].split()[3])])


    # The spacing vector. This will be used in the 'marching_cubes_lewiner'
    # function to define the vertices and faces for 'plot_trisurf' function.
    # This vector is obtained from the sum of the three axis in the cube file as above.
    spacing_vector = axis_1+axis_2+axis_3
    

    # First we read all the isovalues into a 1D list 'isovals'.
    isovals = []

	
    # Starting from the line which the volumetric data starts which is the (natoms+3+2+1+1)th line.
    for i in range(natoms+3+2+1,len(lines)):
        for j in range(0,len(lines[i].split())):
            isovals.append(float(lines[i].split()[j]))
    

    # Define the volumetric numpy array
    wave_fun = np.zeros((nx,ny,nz))
    
    # Setting up the counters to append the isovalues in a 3D numpy array
    c = 0
    c1 = 0
    c2 = 0
	
    for i in range(0,len(isovals)):
        if c2!=nx:
            wave_fun[c2][c1][c] = isovals[i]
            c = c+1
            if c%nz==0:
                c = 0
                c1 = c1+1
                if c1%ny==0:
                    c1 = 0
                    c2 = c2+1
    
    # Now define the x, y, and z 3D arrays to store the grids
    # which is then used for plotting the isosurfaces.
    x_grid = np.zeros((nx,ny,nz))
    y_grid = np.zeros((nx,ny,nz))
    z_grid = np.zeros((nx,ny,nz))
	
    
    # Defining each element of the grid points.
    for i in range(0,nx):
        for j in range(0,ny):
            for k in range(0,nz):
                x_grid[i][j][k] = axis_1[0]*i+axis_2[0]*j+axis_3[0]*k
                y_grid[i][j][k] = axis_1[1]*i+axis_2[1]*j+axis_3[1]*k
                z_grid[i][j][k] = axis_1[2]*i+axis_2[2]*j+axis_3[2]*k



    # For plottin the atoms in the molecule we have to read the
    # coordinates in xyz format which starts from the 7th line.
    coordinates = []
    for i in range(6,natoms+6):
    	coordinates.append(lines[i].split())

    coordinates = np.array(coordinates)


    return coordinates, x_grid, y_grid, z_grid, wave_fun, spacing_vector




def integrate_cube(cube_A, cube_B, grid_volume):
    """
    This function calculates the element-wise multiplication of two numpy arrays 
    and sums their product. Then, it will multiply the sum by 'dv' element to 
    compute the integral of two wavefunction represented as .cube files.
	
    Args:
        cube_A, cube_B (numpy array): The elements of the cube files in a 1D array 
                            obtained from the 'read_cube' function.
        grid_volume (float): The volume of the voxel obtained from grid_volume function.
		
    Returns:
        integral (float): The integration between two wavefunction in the .cube files.
		
    """
    # Compute the element-wise multiplication of the two cube files
    # which were previously obtained in 1D numpy arrays and store 
    # them into another 1D numpy array.
    product = np.multiply(cube_A,cube_B)
    
    # Compute the summation of the above matrix
    summation = product.sum()
    integral = summation*grid_volume
            
    return integral




def form_block_matrix( mat_a, mat_b, mat_c, mat_d ):
    """
    This function gets four numpy arrays and concatenate them into a 
    new matrix in a block format. These matrices should have the same 
    shape on each side which they get concatenated.

         |mat_a   mat_b|
    S =  |             |
         |mat_c   mat_d|

    Args:

        mat_a, mat_b, mat_c, mat_d ( 2D numpy arrays ): The matrices which will form the block matrix

    Returns:

        block_matrix ( numpy 2D array ): The block matrix of the four matrices above.

    """

    # Concatenate the two marix on their row axis
    block_1 = np.concatenate( ( mat_a, mat_b ) )
    block_2 = np.concatenate( (mat_c, mat_d ) )

    # Now concatenate the above concatenated matrices on their 
    # column axis and form the final matrix
    block_matrix = np.concatenate( ( block_1, block_2 ), axis=1 )

    return block_matrix



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

