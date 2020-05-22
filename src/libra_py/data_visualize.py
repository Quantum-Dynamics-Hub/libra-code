#*********************************************************************************                     
#* Copyright (C) 2020 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: data_visualize
   :platform: Unix, Windows
   :synopsis: 
       This module implements various functions for data analysis

.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math
import copy

import numpy as np
from scipy.interpolate import griddata

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn




def plot_map(ax, x_grid, y_grid, z_values, colormap="plasma", resolution=30j):
    """
    This is a function to plot 2D functions

    Args:

        ax ( pyplot instance ): the handler of the plot which we create

        x_grid ( list ): the x grid points, dimension (Nx)

        y_grid ( list ): the y grid points, dimension (Ny)

        z_values ( list of lists ): the values of the function at the grid points, dimension (Nx, Ny)

        colormap ( string ): the type of coloring scheme, 

            Options include: "plasma" (default), "Blues", "viridis", "binary", "hot", etc.

        resolution ( complex, imaginary ): the degree of extra-granulation in the plotting interpolation


    Returns:

        None : just plots the 2D image


    """
    
    
    npts_x = len(x_grid)
    npts_y = len(y_grid)
    
    xmin = x_grid[0]
    xmax = x_grid[npts_x-1]
    
    ymin = y_grid[0]
    ymax = y_grid[npts_y-1]
    
    extent=(xmin, xmax, ymin, ymax)
    
    xs0, ys0, zs0 = [], [], []

    for i in range(npts_x):    
        for j in range(npts_y):
            xs0.append(x_grid[i])
            ys0.append(y_grid[j])
            zs0.append(z_values[i][j])

    #N = 30j
    xs,ys = np.mgrid[extent[0]:extent[1]:resolution, extent[2]:extent[3]:resolution]
    zs = griddata( (xs0, ys0), zs0,  (xs, ys), method="linear")

    #ax.xticks(energy, rotation=30)
    #ax.yticks(energy, rotation=30)    
    
    ax.xticks(rotation=30)
    ax.yticks(rotation=30)
            
    ax.imshow(zs.T, cmap=colormap, extent=extent, interpolation='Lanczos', origin='lower')
    #ax.plot(xs0, ys0, "bo")
    ax.colorbar()

