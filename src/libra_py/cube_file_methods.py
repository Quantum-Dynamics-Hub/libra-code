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


import util.libutil as comn
import libra_py.packages.cp2k.methods as CP2K_methods

def read_cube(filename: str):
    """
    This function reads the wavefunction from a cube file and stores it in
    a 1D numpy array
    Args:
        filename (string): the name of the .cube file to read
    Returns:
        isovalues (numpy.array) : the 1D array of the wavefunctions for all the points on the grid
    """
    
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    # The absolute value is for Gaussian since it might return a negative number for number of atoms
    natoms = abs( int(lines[2].split()[0]) )
    # We skip a few lines in the cube files, that go as follows:
    # 2 lines - comments
    # 1 line  - the number of atoms, etc.
    # 3 lines - the grid spacing and number of grid points in each dimensions

    nstart = natoms+2+1+3        # the index of the first line containing wfc data
    # For Gaussian cube files
    if len( lines[nstart].split() ) < 6:
        nstart += 1

    nlines = len(lines)          # the total number of lines

    isovalues = []
    for i in range(nstart,nlines):
        tmp = lines[i].split()
        ncols = len(tmp)

        for j in range(ncols):
            isovalues.append(float(tmp[j]))
    
    isovalues = np.array(isovalues)
    #data = np.loadtxt(filename,skiprows=n)
    
    return isovalues




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



def plot_cubes( params ):
    """
    This function plots the cubes for selected energy levels using VMD.
    
    Args:
    
        params (dict):
    
            min_band (int): The minimum state number.
                                        
            states_to_be_plotted (list): The list containing the Kohn-Sham orbitals to be plotted by VMD. This list is defined in the submit file.
            
            path_to_tcl_file (str): The path to the tcl file which contains the input for plotting the cubes in VMD.
            
            MO_images_directory (str): The molecular orbitals images directory.
    
            isUKS (int): This parameter is set for spin restricted and unrestricted calculations. When it is
                         set to 1 it means that unrestricted calculations were set in the input file otherwise 
                         it is restricted.
                         
            curr_step (int): The current time step used to save the images of the MOs.
            
            phase_factor_visual (numpy array): The phase correction factor list for each MOs for the current step.
        
    Returns:
    
        None
        
    """
    
   # Critical parameters
    critical_params = [ "min_band", "states_to_be_plotted", "path_to_tcl_file", "MO_images_directory", "curr_step", "phase_factor_visual" ]
    # Default parameters
    default_params = { "isUKS": 0}
    # Check input
    comn.check_input(params, default_params, critical_params)

    # Unpack the Kohn-Sham orbital indicies and the Kohn-Sham orbitals to be plotted. Also unpack the 
    # path to the directory where the molecular orbitals will be plotted
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

    min_band = params["min_band"]
    # We plot the cubes for the previous time step
    cubefile_names_prev = CP2K_methods.cube_file_names_cp2k( params )
    # read the lines of the tcl file
    tcl_file = open(path_to_tcl_file,'r')
    tcl_lines = tcl_file.readlines()
    tcl_file.close()


    if isUKS == 1:
        # The cube file names of the alpha spin is the even indices of the cubefile_names_prev
        alp_cubefile_names_prev = cubefile_names_prev[0::2]
        # The same but for the phase factor of the alpha spin
        phase_factor_alpha      = phase_factor_visual[0::2]
        # The cube file names of the beta spin is the odd indices of the cubefile_names_prev
        bet_cubefile_names_prev = cubefile_names_prev[1::2]
        # The same but for the phase factor of the beta spin
        phase_factor_beta       = phase_factor_visual[1::2]

        for state_to_be_plotted in states_to_be_plotted:
            # Subtracting the min_band to obtain the index of the cube file for the state to be plotted for alpha spin
            alpha_cube_name = alp_cubefile_names_prev[state_to_be_plotted-min_band]
            # Subtracting the min_band to obtain the index of the cube file for the state to be plotted for beta spin
            beta_cube_name = bet_cubefile_names_prev[state_to_be_plotted-min_band]
            # open a new tcl file for alpha cubes
            new_file_alpha = open("vmd_alpha_cube_plot_%d.tcl" % curr_step,'w')
            # open a new tcl file for beta cubes
            new_file_beta = open("vmd_beta_cube_plot_%d.tcl" % curr_step,'w')

            for j in range(len(tcl_lines)):
                if 'mol load cube' in tcl_lines[j]:
                    # Open the cube file in VMD for alpha cubes
                    new_file_alpha.write( 'mol load cube %s\n' % alpha_cube_name )
                    # Open the cube file in VMD for beta cubes
                    new_file_beta.write( 'mol load cube %s\n' % beta_cube_name )
                elif 'render TachyonInternal' in tcl_lines[j]:
                    # Render the images to the MO_images_directory alpha cubes
                    new_file_alpha.write( 'render TachyonInternal %s/%s.tga\n' % ( MO_images_directory, alpha_cube_name.replace('cubefiles/','').replace('.cube','') ) )
                    # Render the images to the MO_images_directory beta cubes
                    new_file_beta.write( 'render TachyonInternal %s/%s.tga\n' % ( MO_images_directory, beta_cube_name.replace('cubefiles/','').replace('.cube','') ) )
                elif 'isosurface' in tcl_lines[j].lower():
                    # Correct the isovalues by multiplying it by the phase factor fo alpha orbitals
                    tmp_elements_alpha = tcl_lines[j].split()
                    tmp_elements_alpha[5] = str( phase_factor_alpha[state_to_be_plotted-min_band] * float( tmp_elements_alpha[5] ) )
                    isosurface_line_alpha = ' '.join(tmp_elements_alpha)
                    new_file_alpha.write( isosurface_line_alpha + '\n' )

                    # Correct the isovalues by multiplying it by the phase factor for beta orbitals
                    tmp_elements_beta = tcl_lines[j].split()
                    tmp_elements_beta[5] = str( phase_factor_beta[state_to_be_plotted-min_band] * float( tmp_elements_beta[5] ) )
                    isosurface_line_beta = ' '.join(tmp_elements_beta)
                    new_file_beta.write( isosurface_line_beta + '\n' )

                else:
                    # The rest of the tcl file lines
                    new_file_alpha.write( tcl_lines[j] )
                    new_file_beta.write( tcl_lines[j] )

            new_file_alpha.close()
            new_file_beta.close()
            # Run the VMD by tcl file
            os.system('vmd < vmd_alpha_cube_plot_%d.tcl' % curr_step)
            os.system('vmd < vmd_beta_cube_plot_%d.tcl' % curr_step)
            #os.system('rm vmd_alpha_cube_plot_%d.tcl' % curr_step)
            #os.system('rm vmd_beta_cube_plot_%d.tcl' % curr_step)

    else:
        # The same as above but with restricted spin calculations. No beta orbitals is considered.
        for state_to_be_plotted in states_to_be_plotted:
            cube_name = cubefile_names_prev[state_to_be_plotted-min_band]

            new_file = open("vmd_cube_plot_%d.tcl" % curr_step,'w')
            for j in range(len(tcl_lines)):
                if 'mol load cube' in tcl_lines[j]:
                    new_file.write( 'mol load cube %s\n' % cube_name )
                elif 'render TachyonInternal' in tcl_lines[j]:
                    new_file.write( 'render TachyonInternal %s/%s.tga\n' % ( MO_images_directory, cube_name.replace('cubefiles/','').replace('.cube','') )  )
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


def plot_cube_v2(params, cube_file_name, phase_factor):
    """
    This function plots the cube files using VMD

    Args:

        params (dictionary):

            vmd_input_template (string): The name of the VMD input template for visualizing the cube files.
     
            states_to_plot (list): The list of states that need to be plot.
     
            plot_phase_corrected (bool): Flag for plotting the molecular orbitals phase-corrected for the running job.
     
            vmd_exe (string): The VMD executable.
     
            tachyon_exe (string): The VMD Tachyon executable for rendering high-quality images.
     
            x_pixels (integer): Number of pixels in the X direction of the image.
     
            y_pixels (integer): Number of pixels in the Y direction of the image.
     
            remove_cube (bool): Flag for removing the cube files after plotting the molecular orbitals.

        cube_file_name (string): The name of the cube file

        phase_factor (integer): The phase factor for plotting the phase-corrected molecular orbital

    Returns: 

        None
    """

    print(F'Plotting cube file {cube_file_name}...')
    file = open(params['vmd_input_template'], 'r')
    tcl_lines = file.readlines()
    file.close()
    vmd_exe = params['vmd_exe']
    tachyon_exe = params['tachyon_exe']
    x_pixels = params['x_pixels']
    y_pixels = params['y_pixels']
    image_format = params['image_format']
    together_mode = params['together_mode']
    if together_mode:
        new_tcl_name = 'vmd_tmode.tcl'
        file = open(new_tcl_name, 'w')
        cube_file_names = []
        for state in params['states_to_plot']:
            cube_file_names.append(cube_file_name.split('WFN')[0]+F'WFN_{str(state).zfill(5)}_1-1_0.cube')
        print(cube_file_names)
        state_counter = 0
        for i in range(len(tcl_lines)):
            if 'load cube' in tcl_lines[i]:
                tmp_name = cube_file_names[state_counter]
                file.write(F'mol load cube {tmp_name}\n')
                print(tmp_name)
                state_counter += 1
            elif 'render' in tcl_lines[i]:
                tmp_name = cube_file_name.split('WFN')[0]
                file.write(F'render Tachyon {tmp_name} "{tachyon_exe} -aasamples 12 %s -format {image_format.upper()} -res {x_pixels} {y_pixels} -o %s.{image_format.lower()}"\n')
            else:
                file.write(tcl_lines[i])
    
        file.close()
    else:
        state_name = cube_file_name.replace('.cube','')
        new_tcl_name = state_name+'.tcl'
        file = open(new_tcl_name, 'w')
    
        for i in range(len(tcl_lines)):
            if 'load cube' in tcl_lines[i]:
                file.write(F'mol load cube {cube_file_name}\n')
            elif 'Isosurface' in tcl_lines[i]:
                tmp = tcl_lines[i].split()
                for k in range(len(tmp)):
                    if tmp[k]=='Isosurface':
                        break
                tmp[k+1] = str(float(tmp[k+1]) * phase_factor)
                tcl_lines[i] = ' '.join(tmp) + '\n'
                file.write(tcl_lines[i])
            elif 'render' in tcl_lines[i]:
                file.write(F'render Tachyon {state_name} "{tachyon_exe} -aasamples 12 %s -format {image_format.upper()} -res {x_pixels} {y_pixels} -o %s.{image_format.lower()}"\n')
            else:
                file.write(tcl_lines[i])
    
        file.close()

    os.system(F'{vmd_exe} < {new_tcl_name}')
    if params['remove_cube']:
        if together_mode:
            for name in cube_file_names:
                os.system(F'rm {name}')
            os.system(F'rm {tmp_name}')
        else:
            os.system(F'rm {cube_file_name}')
            os.system(F'rm {state_name}')
    




