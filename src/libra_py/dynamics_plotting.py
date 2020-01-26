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
.. moduleauthor:: Alexey V. Akimov

  List of functions:
  
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


def plot_dyn(plot_params):
    """
    This function is meant to plot the results stored in the hdf files generated by the dynamics runs

    Args:

        prefix ( string ): the name of the directory containing the input HDF5 file
            This directory will also be used to output the generated picture files [ default : "out"]
        filename ( string ): name of the HDF5 file to read [ default: "data.hdf"]
        output_level ( int ): the level of info contained in the HDF5 file [ default : 3]
        which_trajectories ( list of ints ) : indices of the trajectories to print [ default: [0] ]
        which_adi_states ( list of ints ) : indices of the adiabatic states to print [ default: [0] ]
        which_dia_states ( list of ints ) : indices of the diabatic states to print [ default: [0] ]
        colors ( dictionary ): the definition of the colors to use
        clrs_index ( list of strings ) : defines the mapping of the colors on integers and vice versa 

    See in the ```dynamics_hdf5.py``` file for the data sets available                    

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
    default_params = {  "prefix":"out", "filename":"data.hdf", "output_level":3,
                        "which_trajectories":[0], "which_dofs":[0],
                        "which_adi_states":[0], "which_dia_states":[0],
                        "colors":colors, "clrs_index":clrs_index
                     }
    comn.check_input(plot_params, default_params, critical_params)
        

    filename = plot_params["filename"]
    prefix = plot_params["prefix"]
    output_level = plot_params["output_level"]
    which_trajectories = plot_params["which_trajectories"]
    which_dofs = plot_params["which_dofs"]
    which_adi_states = plot_params["which_adi_states"]
    which_dia_states = plot_params["which_dia_states"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]

    out_prefix = prefix

 
    with h5py.File(F"{prefix}/{filename}", 'r') as f:

        _figsize = (24,24)

        if output_level>=1:
            
            #========= Energies and their fluctuations =========

            plt.figure(num=None, figsize=(24, 12), dpi=300, frameon=False)        
            
            plt.subplot(1,2,1)            
            plt.title('t-Energies')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Energy, a.u.')        
            
            plt.plot(f["time/data"][:], f["Ekin_ave/data"][:], label="Kinetic energy", linewidth=2, color = colors["11"])                 
            plt.plot(f["time/data"][:], f["Epot_ave/data"][:], label="Potential energy", linewidth=2, color = colors["21"])                 
            plt.plot(f["time/data"][:], f["Etot_ave/data"][:], label="Total energy", linewidth=2, color = colors["31"])                 
            
            
            plt.subplot(1,2,2)            
            plt.title('t-dEnergies')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Energy fluctuation, a.u.')        
            
            plt.plot(f["time/data"][:], f["dEkin_ave/data"][:], label="Kinetic energy", linewidth=2, color = colors["11"])                 
            plt.plot(f["time/data"][:], f["dEpot_ave/data"][:], label="Potential energy", linewidth=2, color = colors["21"])                 
            plt.plot(f["time/data"][:], f["dEtot_ave/data"][:], label="Total energy", linewidth=2, color = colors["31"])                             
                
                
            plt.savefig(F"{out_prefix}/t-en.png", dpi=300)
            plt.show()
            

        if output_level>=2:
            
            #========= State indices =========
            
            #ntraj = f["states"].attrs["dim"][0]            
            ntraj = f["time/data"].shape[0]

            plt.figure(num=None, figsize=(24, 24), dpi=300, frameon=False)        
                        
            plt.title('State indices')
            plt.xlabel('Time, a.u.')
            plt.ylabel('State index')        
            
            for tr in range(ntraj):
                if tr in which_trajectories:
                    plt.plot(f["time/data"][:], f["states/data"][:, tr], label="", linewidth=2, color = colors["11"])                             
                            
            plt.savefig(F"{out_prefix}/t-state_indices.png", dpi=300)
            plt.show()
            
            
            
        if output_level>=3:
          
            #========= SH, adi SE and dia SE populations =========
            nadi = f["D_adi/data"].shape[1] #attrs['dim'][1]
            ndia = f["D_dia/data"].shape[1] #attrs['dim'][1]

            plt.figure(num=None, figsize=(24, 8), dpi=300, frameon=False)        

            #================ SH populations =============
            plt.subplot(1,3,1)            
            plt.title('SH populations')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Population, a.u.')        
                        
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    indx = indx + 1
                    plt.plot(f["time/data"][:], f["SH_pop/data"][:, istate, 0], label=F"state {istate}", linewidth=2, color = colors[clrs_index[indx] ])                 

            #================ adi SE populations =============
            plt.subplot(1,3,2)            
            plt.title('adi SE populations')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Population, a.u.')        
            
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    indx = indx + 1
                    plt.plot(f["time/data"][:], f["D_adi/data"][:, istate, istate], label=F"state {istate}", linewidth=2, color = colors[clrs_index[indx] ])                 
                    
                    
            #================ dia SE populations =============
            plt.subplot(1,3,3)            
            plt.title('dia SE populations')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Population, a.u.')        
            
            indx = -1
            for istate in range(ndia):
                if istate in which_dia_states:
                    indx = indx + 1
                    plt.plot(f["time/data"][:], f["D_dia/data"][:, istate, istate], label=F"state {istate}", linewidth=2, color = colors[clrs_index[indx] ])                 
                    
                    
            plt.savefig(F"{out_prefix}/t-pops.png", dpi=300)
            plt.show()            
        
        
        
        if output_level>=3:
            
            #========= Coordinates =========
            ntraj = f["q/data"].shape[1] #attrs['dim'][1]
            ndof = f["q/data"].shape[2] #attrs['dim'][2]
        

            plt.figure(num=None, figsize=(24, 12), dpi=300, frameon=False)        
            plt.subplot(1,2,1)            
            plt.title('t-q')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Coordiante, a.u.')        
            
            indx = -1
            for idof in range(ndof):
                if idof in which_dofs:
                    indx = indx + 1
                    for tr in range(ntraj):
                        if tr in which_trajectories:
                            plt.plot(f["time/data"][:], f["q/data"][:, tr, idof], label="", linewidth=2, color = colors[ clrs_index[indx] ])                 
                            
                        
            #========= Phase-space =========            
            plt.subplot(1,2,2)            
            plt.title('q-p')        
            plt.xlabel('Coordiante, a.u.')        
            plt.ylabel('Momenta, a.u.')
            
            indx = -1
            for idof in range(ndof):
                if idof in which_dofs:
                    indx = indx + 1
                    for tr in range(ntraj):
                        if tr in which_trajectories:                                    
                            plt.plot(f["q/data"][:, tr, idof], f["p/data"][:, tr, idof], label="", linewidth=2, color = colors[ clrs_index[indx] ])                 
                
            plt.savefig(F"{out_prefix}/t-q-p.png", dpi=300)
            plt.show()
            
            
        if output_level>=4:            

            ntraj = f["hvib_adi/data"].shape[1] #attrs['dim'][1]
            nadi = f["hvib_adi/data"].shape[2] #attrs['dim'][2]
            ndia = f["hvib_dia/data"].shape[2] # attrs['dim'][2]

            #============== Adiabatic energies =============
            plt.figure(num=None, figsize=(24, 12), dpi=300, frameon=False)        
            plt.subplot(1,2,1)            
            plt.title('t-H_adi')        
            plt.xlabel('Time, a.u.')        
            plt.ylabel('Energy, a.u.')
            
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    indx = indx + 1
                    for tr in range(ntraj):
                        if tr in which_trajectories:                                    
                            plt.plot(f["time/data"][:], f["hvib_adi/data"][:, tr, istate, istate], label="", linewidth=2, color = colors[ clrs_index[indx] ])                 

            #============== Diabatic energies =============                            
            plt.subplot(1,2,2)            
            plt.title('t-H_dia')        
            plt.xlabel('Time, a.u.')        
            plt.ylabel('Energy, a.u.')
            
            indx = -1
            for istate in range(ndia):
                if istate in which_dia_states:
                    indx = indx + 1
                    for tr in range(ntraj):
                        if tr in which_trajectories:                                    
                            plt.plot(f["time/data"][:], f["hvib_dia/data"][:, tr, istate, istate], label="", linewidth=2, color = colors[ clrs_index[indx] ])
                            
            plt.savefig(F"{out_prefix}/t-hvib.png", dpi=300)
            plt.show()

            
        if output_level>=4:            

            ntraj = f["St/data"].shape[1] #attrs['dim'][1]
            nadi = f["St/data"].shape[2] #attrs['dim'][2]            

            #============== St diagonal matrix elements =============
            plt.figure(num=None, figsize=(24, 12), dpi=300, frameon=False)        
            plt.subplot(1,2,1)            
            plt.title('t-St')        
            plt.xlabel('Time, a.u.')        
            plt.ylabel('Time overlap')
            
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    indx = indx + 1
                    for tr in range(ntraj):
                        if tr in which_trajectories:                                    
                            plt.plot(f["time/data"][:], f["St/data"][:, tr, istate, istate], label="", linewidth=2, color = colors[ clrs_index[indx] ])                 

            #============== Projectors =============                            
            plt.subplot(1,2,2)            
            plt.title('t-projectors')        
            plt.xlabel('Time, a.u.')        
            plt.ylabel('Projectors')
            
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    indx = indx + 1
                    for tr in range(ntraj):
                        if tr in which_trajectories:                                    
                            plt.plot(f["time/data"][:], f["projector/data"][:, tr, istate, istate], label="", linewidth=2, color = colors[ clrs_index[indx] ])
                            
            plt.savefig(F"{out_prefix}/St-projector.png", dpi=300)
            plt.show()
                            

        if output_level>=4:            

            ntraj = f["basis_transform/data"].shape[1] #attrs['dim'][1]
            ndia  = f["basis_transform/data"].shape[2] #attrs['dim'][2]            
            nadi  = f["basis_transform/data"].shape[3] #attrs['dim'][3]            

            #============== dia-to-adia transformation matrix =============
            plt.figure(num=None, figsize=(24, 24), dpi=300, frameon=False)        
            plt.subplot(1,1,1)            
            plt.title('t-Basis_transform')        
            plt.xlabel('Time, a.u.')        
            plt.ylabel('Basis transform')
            
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    for istate2 in range(ndia):
                        if istate2 in which_dia_states:
                            
                            indx = indx + 1
                            for tr in range(ntraj):
                                if tr in which_trajectories:                                    
                                    plt.plot(f["time/data"][:], f["basis_transform/data"][:, tr, istate2, istate], label="", linewidth=2, color = colors[ clrs_index[indx] ])                 
            
            plt.savefig(F"{out_prefix}/basis_transform.png", dpi=300)
            plt.show()
                                            
                
        plt.close()
        



def plot_dyn_old(res):
    """
    This is a function to plot dynamical results for 2-state system
    Yeah, this is not a general-purpose function, but is a useful exmaple
  
    Format: 
       0      1      2       3        4        5           6          7        8            9       10        11            12        13       14         15             16          17
    #obs_T, obs_q, obs_p, obs_Ekin, obs_Epot, obs_Etot, obs_dEkin, obs_dEpot, obs_dEtot, obs_Cadi, obs_Cdia, obs_dm_adi, obs_dm_dia, obs_pop  obs_states obs_hvib_adi obs_hvib_dia  obs_St
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



    obs_T = res[0]
    obs_Ekin = res[3]
    obs_Epot = res[4]
    obs_Etot = res[5]
    obs_dEkin = res[6]
    obs_dEpot = res[7]
    obs_dEtot = res[8]
    obs_dm_adi00 = data_conv.unpack1(res[11], 0, 0, 0)
    obs_dm_adi11 = data_conv.unpack1(res[11], 1, 1, 0)
    obs_dm_dia00 = data_conv.unpack1(res[12], 0, 0, 0)
    obs_dm_dia11 = data_conv.unpack1(res[12], 1, 1, 0)
    obs_pop00 = data_conv.unpack1(res[13], 0, 0, 2)
    obs_pop11 = data_conv.unpack1(res[13], 1, 0, 2)
    obs_q = data_conv.unpack1(res[1], 0, 0, 2)
    obs_p = data_conv.unpack1(res[2], 0, 0, 2)
    
    ndof = res[1][0].num_of_rows
    ntraj = res[1][0].num_of_cols
    

    plt.figure(1, figsize=(24, 24)) # dpi=300, frameon=False)

    plt.subplot(3,2,1)
    plt.title('q-t')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Coordiante, a.u.')
    
    for tr in range(ntraj):
        obs_q = data_conv.unpack1(res[1], 0, tr, 2)
        plt.plot(obs_T, obs_q, label="", linewidth=2, color = colors["11"])   # label='q vs. t'
    plt.legend()

    plt.subplot(3,2,2)
    plt.title('Phase portrait')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Momentum, a.u.')    
    for tr in range(ntraj):
        obs_q = data_conv.unpack1(res[1], 0, tr, 2)
        obs_p = data_conv.unpack1(res[2], 0, tr, 2)
        plt.plot(obs_q, obs_p, linewidth=2, color = colors["11"])  # label='p vs. q'
    plt.legend()

    plt.subplot(3,2,3)
    plt.title('Adiabatic populations')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population of state 0')
    plt.plot(obs_T, obs_pop00, "--o", label='SH P(0)', linewidth=2, color = colors["11"]) 
    plt.plot(obs_T, obs_dm_adi00, "--o", label='SE P(0)', linewidth=2, color = colors["21"]) 
    plt.plot(obs_T, obs_pop11, "--", label='SH P(1)', linewidth=2, color = colors["11"]) 
    plt.plot(obs_T, obs_dm_adi11, "--", label='SE P(1)', linewidth=2, color = colors["21"]) 
    plt.legend()

    plt.subplot(3,2,4)
    plt.title('Diabatic populations')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population of state 0')
    plt.plot(obs_T, obs_dm_dia00, label='SE P(0)', linewidth=2, color = colors["21"]) 
    plt.plot(obs_T, obs_dm_dia11, "--", label='SE P(1)', linewidth=2, color = colors["21"]) 
    plt.legend()

    plt.subplot(3,2,5)
    plt.title('Energies')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(obs_T, obs_Ekin, label='Ekin', linewidth=2, color = colors["11"]) 
    plt.plot(obs_T, obs_Epot, label='Epot', linewidth=2, color = colors["21"]) 
    plt.plot(obs_T, obs_Etot, label='Etot', linewidth=2, color = colors["31"]) 
    plt.legend()

    plt.subplot(3,2,6)
    plt.title('Energy fluctuations')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(obs_T, obs_dEkin, label='dEkin', linewidth=2, color = colors["11"]) 
    plt.plot(obs_T, obs_dEpot, label='dEpot', linewidth=2, color = colors["21"]) 
    plt.plot(obs_T, obs_dEtot, label='dEtot', linewidth=2, color = colors["31"]) 
    plt.legend()

    plt.show()
    plt.close()

