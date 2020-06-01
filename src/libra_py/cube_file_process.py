#*********************************************************************************
#* Copyright (C) 2020 Mohammad Shakiba, Brendan Smith, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: cube_file_process
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing with cube files
.. moduleauthors:: 
       Mohammad Shakiba, Brendan Smith, Alexey V. Akimov 
  
"""


import os
import sys
import math
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import numpy as np




def read_cube(filename: str):
    """
    This function reads the wavefunction from a cube file and stores it in
    a 1D numpy array

    Args:

        filename ( string ): the name of the .cube file to read

    Returns:

        numpy.array : the 1D array of the wavefunctions for all the points on the grid
    

    Note: 

        Previously, it was called `read_1`

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
		
        Note:
        
	    It uses the formula of calculating the parallelpiped volume:
            https://en.wikipedia.org/wiki/Parallelepiped
        	
    """
	
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    
    # We use the 3rd, 4th and 5th row in the lines to obtain
    # the vectors of the parallelpiped into a numpy array (V).
    V = np.array([[float(lines[3].split()[1]),float(lines[3].split()[2]),float(lines[3].split()[3])],\
                 [float(lines[4].split()[1]),float(lines[4].split()[2]),float(lines[4].split()[3])],\
                 [float(lines[5].split()[1]),float(lines[5].split()[2]),float(lines[5].split()[3])]\
                 ])
    # Then we calculate the determinant of the a to obtain the volume.
    dv = np.absolute(np.linalg.det(V))
	
    
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
	
        Coordinates (2D numpy array): The coordinates of the molecule structure in
                                      the same format shown in the .cube file.
		
           x, y, z (3D numpy array): Containing the grid points for each of the X-,
                                     Y-, and Z- axis.
		
           v (3D numpy array): The volumetric data in a 3D numpy array format.
		
           spacing (numpy 1D array): The spacing vector used for plotting the isosurfaces.

    """

    f = open(filename,'r')
    lines = f.readlines()
    f.close()
	
	# The number of atoms in the 3rd line
    natoms = int(lines[2].split()[0])
	
	# The number of voxels defined for each axis
    vx = int(lines[3].split()[0])
    vy = int(lines[4].split()[0])
    vz = int(lines[5].split()[0])
	
	# The spacing vector. This will be used in the 'marching_cubes_lewiner'
	# function to define the vertices and faces for 'plot_trisurf' function.
    spacing = []
    spacing.append([float(lines[3].split()[1])+float(lines[4].split()[1])+float(lines[5].split()[1]),\
               float(lines[3].split()[2])+float(lines[4].split()[2])+float(lines[5].split()[2]),\
               float(lines[3].split()[3])+float(lines[4].split()[3])+float(lines[5].split()[3])])
    
    # Bohr to Angstrom : *0.529177
	# Here we use the same unit as is used in the .cube file structure which is Bohr.
	# The three vectors below are the same vectors which are present in lines 4th to 6th
	# which are used to create the 3D numpy array of the grid points (x, y, z) which are
	# used to plot the isosurfaces.
    vecx = np.array([float(lines[3].split()[1]),float(lines[3].split()[2]),float(lines[3].split()[3])])#*0.529177
    vecy = np.array([float(lines[4].split()[1]),float(lines[4].split()[2]),float(lines[4].split()[3])])#*0.529177
    vecz = np.array([float(lines[5].split()[1]),float(lines[5].split()[2]),float(lines[5].split()[3])])#*0.529177

    # First we read all the isovalues in a 1D array
    s = []

    for i in range(natoms+3+2+1,len(lines)):
        for j in range(0,len(lines[i].split())):
            s.append(float(lines[i].split()[j]))

    v = np.zeros((vx,vy,vz)) # Define the volumetric numpy matrix
    
	# Setting up the counters
    c = 0
    c1 = 0
    c2 = 0
	
    for i in range(0,len(s)):
        if c2!=vx:
            v[c2][c1][c] = s[i]
            c = c+1
            if c%vz==0:
                c = 0
                c1 = c1+1
                if c1%vy==0:
                    c1 = 0
                    c2 = c2+1
    
	# Now define the x, y, and z 3D arrays to store the grids
	# which is then used for plotting the isosurfaces.
    x = np.zeros((vx,vy,vz))
    y = np.zeros((vx,vy,vz))
    z = np.zeros((vx,vy,vz))
	
	# Defining each element of the grid points.
    for i in range(0,vx):
        for j in range(0,vy):
            for k in range(0,vz):
                x[i][j][k] = vecx[0]*i+vecy[0]*j+vecz[0]*k
                y[i][j][k] = vecx[1]*i+vecy[1]*j+vecz[1]*k
                z[i][j][k] = vecx[2]*i+vecy[2]*j+vecz[2]*k

    # For plottin the atoms in the molecule we have to read the
	# coordinates in xyz format which starts from the 7th line.
    coordinates = []
    for i in range(6,natoms+6):
    	coordinates.append(lines[i].split())

    coordinates = np.array(coordinates)
                
    return coordinates, x, y, z, v, spacing




def cube_ovlp_arb(A, B, dv):
    """
    This function calculates the element-wise multiplication of two numpy arrays 
    and sums their product. Then, it will multiply the sum by 'dv' element to 
    compute the integral of two wavefunction represented as .cube files.
	
    Args:
        A, B (numpy array): The elements of the cube files in a 1D array 
                            obtained from the 'read_cube' function.
        dv (float): The volume of the voxel obtained from grid_volume function.
		
    Returns:
        S*dv (float): The integration between two wavefunction in the .cube files.
		
    """

    C = np.multiply(A,B)
    S = C.sum()
            
    return S*dv


