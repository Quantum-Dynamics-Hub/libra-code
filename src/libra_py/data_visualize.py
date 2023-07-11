#*********************************************************************************                     
#* Copyright (C) 2020 Brendan Smith, Alexey V. Akimov                                                   
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

.. moduleauthor:: Alexey V. Akimov, Brendan Smith

"""

import os
import sys
import math
import copy
import scipy

import numpy as np
from scipy.interpolate import griddata

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn




colors = {}

colors.update({"11": "#8b1a0e"})  # red       
colors.update({"12": "#FF4500"})  # orangered 
colors.update({"13": "#B22222"})  # firebrick 
colors.update({"14": "#DC143C"})  # crimson   

colors.update({"21": "#5e9c36"})  # green
colors.update({"22": "#006400"})  # darkgreen  
colors.update({"23": "#228B22"})  # forestgreen
colors.update({"24": "#808000"})  # olive      

colors.update({"31": "#8A2BE2"})  # blueviolet
colors.update({"32": "#00008B"})  # darkblue  

colors.update({"41": "#2F4F4F"})  # darkslategray

clrs_index = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]



def plot_map(ax, x_grid, y_grid, z_values, colormap="plasma", resolution=30j, savefig=0, figure_name=None):
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

        savefig ( int ):  0 - don't save the figure as a file, 1 - do save it

        figure_name ( string ): the name of the file to where the figure is to be saved (only is used if savefig == 1)


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

    if savefig == 1:
        ax.savefig(figure_name)





def plot_map_nparray( plt, data, colormap='hot', interpolation_scheme='nearest', fig_width=6.42, fig_height=2.14, num_subplots=1, \
                      titles=[""], title_size=10, xticks_size=8, yticks_size=8, colorbar_label='fs',
                      savefig=0, showfig=1, figure_names=['plotted_map_nparray.png'] ):

    """
    This is a function to plot a heat map using a 2d numpy array

    Args:

        plt ( pyplot instance ): the handler of the plot which we create

        data ( list of  2d numpy arrays ): Contains the data to be plotted using a heat map. Each list element can be its own subplot
        
        colormap ( string ): The type of coloring scheme, 

            Options include: "plasma", "Blues", "viridis", "binary", "hot" [default] , etc.
        
        interpolation_scheme ( string ): Describes how to terpolate the presented data. 
                                         Nearest for no interpolation [ defaut: "nearest" ]
        
        fig_width ( float ): width of the final figure. May contain multiple subplots [ default: 6.42 ]
        
        fig_height ( float ): height of the final figure. May contain multiple subplots [ default: 2.14 ]

        num_subplots ( int ): number of subplots to have in the total figure [ default: 1 ]

        titles ( list of strings ): the title for each subplot [ default: [""] ]

        title_size ( int ): fontsize of the title [ default: 10 ]
        
        xticks_size ( int ): fontsize of the title [ default: 8 ]
        
        yticks_size ( int ): fontsize of the title [ default: 8 ]
        
        colorbar_label ( string ): Name of the label of the colorbar [ default: "fs" ]

        savefig ( int ): flag telling whether to save figure or not; 0 - no, 1 - yes [ default: 0 ]

        showfig ( int ): flat telling whether to show figure or not; 0 - no, 1 - yes [ default: 1 ]

        figure_names ( list of strings ): names of the files for all subplots [ default: "plotted_map_nparray.png" ] 

    Returns:

        None : just plots the 2D image(s)
    """
    
    if len(titles) != num_subplots:
        print("\nYou do not have a title for every subplot. Assigning the title for subplot[0] to all subplots")
        for i in range(num_subplots-1):
            titles.append( titles[0] )

    if len(data) != num_subplots:
        print("\nYou do not have a data for every subplot. Assigning the data for subplot[0] to all subplots")
        for i in range(num_subplots-1):
            data.append( data[0] )
    
    plt.figure(num=None, figsize=(fig_width, fig_height), dpi=300)
    for sub_plot in range( num_subplots ):

        plt.subplot(1, num_subplots, sub_plot+1)
        plt.title(titles[sub_plot], fontsize=title_size)
        plt.xticks(fontsize=xticks_size)
        plt.yticks(fontsize=yticks_size)
        plt.imshow( data[sub_plot], cmap=colormap, interpolation=interpolation_scheme )
        cb = plt.colorbar(label=colorbar_label)
        cb.ax.tick_params(labelsize=8)
        plt.tight_layout()

        if savefig == 1:
            plt.savefig(figure_names[sub_plot])

        if showfig == 1:
            plt.show()


def plot_nparrays( plt, xdata, ydata, fig_width=6.42, fig_height=2.14, num_subplots=1, colors=["black"], \
                   titles=[""], labels=[ [""] ], title_size=10, xticks_size=8, yticks_size=8, legend_fontsize=8, \
                   n_superimposing=1, savefig=0, figure_names=['plotted_nparray.png']):

    """
    This function takes in n nparrays and plots them using the features described below


    Args:

        plt ( pyplot instance ): the handler of the plot which we create

        xdata ( list of numpy arrays ): Contains the data to be plotted. Each list element can be its own subplot

        ydata ( list of numpy arrays ): Contains the data to be plotted. Each list element can be its own subplot

        fig_width ( float ): width of the final figure. May contain multiple subplots 
        
        fig_height ( float ): height of the final figure. May contain multiple subplots         

        num_subplots ( int ): number of subplots to have in the total figure        

        titles ( list of strings ): the title for each subplot     

        title_size ( int ): fontsize of the title
        
        xticks_size ( int ): fontsize of the title
        
        yticks_size ( int ): fontsize of the title
        
        n_superimposing ( int ): number of data to super impose on one plot

    Returns:

        None : just plots the 2D image(s)
    """

    if len(titles) != num_subplots:
        print("\nYou do not have a title for every subplot. Assigning titles[0] to all subplots")
        for i in range(num_subplots-1):
            titles.append( titles[0] )

    if len(labels) != num_subplots:
        print("\nYou do not have a labels for every subplot. Assigning labels[0] to all subplots")
        for i in range(num_subplots-1):
            labels.append( labels[0] )

    if len(xdata) != num_subplots:
        print("\nYou do not have a xdata for every subplot. Assigning xdata[0] to all subplots")
        for i in range(num_subplots-1):
            xdata.append( xdata[0] )

    if len(ydata) != num_subplots:
        print("\nYou do not have a ydata for every subplot. Assigning ydata[0] to all subplots")
        for i in range(num_subplots-1):
            ydata.append( ydata[0] )

    if len(figure_names) != num_subplots:
        print("\nYou do not have a figure_name for every subplot. Assigning a figure_name to all subplots")
        for i in range(num_subplots-1):
            figure_names.append( "tmp_"+str(i)+".png" )

    plt.figure(num=None, figsize=(fig_width, fig_height), dpi=300)
    for sub_plot in range( num_subplots ):

        plt.subplot(1, num_subplots, sub_plot+1)
        plt.title(titles[sub_plot], fontsize=title_size)
        plt.xticks(fontsize=xticks_size)
        plt.yticks(fontsize=yticks_size)

        for superimpose in range( n_superimposing ):

            plt.plot( xdata[sub_plot][0], ydata[sub_plot][superimpose], label=labels[sub_plot][superimpose] )

        plt.legend(fontsize=legend_fontsize)
        plt.tight_layout()

        if savefig == 1:
            plt.savefig(figure_names[sub_plot])
 
