#*********************************************************************************                     
#* Copyright (C) 2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: dynamics_plotting
   :platform: Unix, Windows
   :synopsis: This module implements the functions to read the HDF5 files with the results of the dynamical calculations
       List of functions:
           * plot_surfaces(_compute_model, _param_sets, states_of_interest, xmin, xmax, dx, plot_params)
           * plot_dyn(plot_params)
           * plot_dyn_old(res)

.. moduleauthor:: Alexey V. Akimov
  
"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019 Alexey V. Akimov"
__credits__ = ["Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"


import h5py
import numpy as np

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn

import matplotlib.pyplot as plt   # plots
#from matplotlib.mlab import griddata


def plot_surfaces(_compute_model, _param_sets, states_of_interest, xmin, xmax, dx, plot_params):
    """
    Args:
        _compute_model ( PyObject ): the function that returns the class with Hamiltonian properties
        _param_sets ( list of lists of dictionaries): parameters of the models
        states_of_interest ( list of ints ): indices of the states we want to plot 
        
    """


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



    # Parameters and dimensions
    critical_params = [  ] 
    default_params = {  "colors":colors, "clrs_index":clrs_index,
                        "xlim":[-7.5, 15], "ylim":[-0.005, 0.025]
                     }
    comn.check_input(plot_params, default_params, critical_params)
        
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]

    
    X = []
    nsteps = int((xmax - xmin) / dx) + 1

    for i in range(nsteps):
        X.append(xmin + i * dx)
    
    

    plt.rc('axes', titlesize=38)      # fontsize of the axes title
    plt.rc('axes', labelsize=38)      # fontsize of the x and y labels
    plt.rc('legend', fontsize=36)     # legend fontsize
    plt.rc('xtick', labelsize=28)     # fontsize of the tick labels
    plt.rc('ytick', labelsize=28)     # fontsize of the tick labels

    plt.rc('figure.subplot', left=0.2)
    plt.rc('figure.subplot', right=0.95)
    plt.rc('figure.subplot', bottom=0.13)
    plt.rc('figure.subplot', top=0.88)
    

    sz = len(_param_sets)
    for iset in range(sz):
        
        n = len(_param_sets[iset]["E_n"])
        
        ham = nHamiltonian(n, n, 1) # ndia, nadi, nnucl
        ham.init_all(2)

        
        hdia, hadi  = [], []  
        uij = []              # projecitions of the MOs onto elementary basis
    
        for k1 in range(n):
            hadi.append([])
            hdia.append([])
            uij_k1 = []
            for k2 in range(n):
                uij_k1.append([])
            uij.append(uij_k1)
                
    
        for i in range(nsteps):
                 
            q = MATRIX(1,1); q.set(0, 0, X[i])        
        
            # Diabatic properties
            ham.compute_diabatic(_compute_model, q, _param_sets[iset])
            
            # Adiabatic properties
            ham.compute_adiabatic(1);       
        
            U = ham.get_basis_transform()
            #P = U * U.H()  # population matrix
        
            for k1 in range(n):
                hadi[k1].append(ham.get_ham_adi().get(k1,k1).real)
                hdia[k1].append(ham.get_ham_dia().get(k1,k1).real)
            
                for k2 in range(n):
                    uij[k1][k2].append(U.get(k1,k2).real**2 + U.get(k1,k2).imag**2)
                    

                    
        plt.figure(2*iset, figsize=(36, 18)) # dpi=300, frameon=False)
                
        plt.subplot(1, 2, 1)    
        plt.ylim(ylim[0], ylim[1])
        plt.xlim(xlim[0], xlim[1])
    
        #plt.title('Params set %i: Ham_dia' % (iset) )
        plt.xlabel('Coordinate, a.u.')
        plt.ylabel('Energy, a.u.')
        for k1 in range(n):
            plt.plot(X, hdia[k1], label='$H_{%i%i}$' % (k1,k1), linewidth=7, color = colors[clrs_index[k1]])     
        plt.legend()    
    
        plt.subplot(1, 2, 2)
        plt.ylim(ylim[0], ylim[1])
        plt.xlim(xlim[0], xlim[1])
        #plt.title('Params set %i: Ham_adi' % (iset))
        plt.xlabel('Coordinate, a.u.')
        plt.ylabel('Energy, a.u.')
        for k1 in range(n):
            plt.plot(X, hadi[k1], label='$E_{%i}$' % (k1), linewidth=7, color = colors[clrs_index[k1]])     
        plt.legend()    
            
            
          
        plt.figure(2*iset+1, figsize=(36, 18)) # dpi=300, frameon=False)            
        sz1 = len(states_of_interest)
    
        for k2 in states_of_interest:
            indx = states_of_interest.index(k2)
            plt.subplot(1, sz1, 1+indx)
            #plt.title('Params set %i: Adi state %i' % (iset, k2) )
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('Projection')
    
            for k1 in range(n):            
                plt.plot(X, uij[k1][k2], label='$<%i|%i>$' % (k1, k2), linewidth=7, color = colors[clrs_index[k1]])         
            plt.legend()               
    

        plt.show()
        plt.close()          

