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



