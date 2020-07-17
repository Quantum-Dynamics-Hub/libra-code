#*********************************************************************************                     
#* Copyright (C) 2019-2020 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: plotters
   :platform: Unix, Windows
   :synopsis: This module implements functions for visualizing results of the exact dynamics calculations

       List of functions:
           * exact_init_hdf5(saver, hdf5_output_level, _nsteps, _ndof, _nstates, _ngrid)
           * exact_init_custom_hdf5(saver, _nsteps, _ncustom_pops, _nstates)
           * save_data_hdf5(step, wfc, saver, params)
           * save_data_mem(step, wfc, saver, params)


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


import os
import sys
import math
import copy

import h5py
import numpy as np
import matplotlib.pyplot as plt   # plots

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
import libra_py.units as units
import libra_py.data_outs as data_outs
import libra_py.data_conv as data_conv
import libra_py.data_savers as data_savers




def plot_mem(res, _params, model_params, plot_params):
    """
    Args:
        res ( mem_saver ): the object with all the suitable parameters
        _params ( dict ): simulation control parameters
        model_params ( dict ): the parameters of the model we compute
        
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
    default_params = {  "colors":colors, "clrs_index":clrs_index }
    comn.check_input(plot_params, default_params, critical_params)
        
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]



    
    
    params = dict(_params)
    
    nsteps = params["nsteps"]    
    nstates = len(model_params["E_n"])
    prefix = params["prefix"]
    properties_to_save = params["properties_to_save"]
 
    
    t = None
    if "time" in properties_to_save:
        #t = res.data["time"]
        t = list(res.np_data["time"][:])
    
    #=============== Populations ======================
    
    plt.figure(1, figsize=(36, 12)) # dpi=300, frameon=False)    
    plt.subplot(1, 2, 1)
    plt.title('Adiabatic population dynamics' )
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    
    if "pop_adi" in properties_to_save and t != None:
        for i in range(nstates):        
            #Pi = data_conv.unpack1(res.data["pop_adi"], i, 0, 0)                
            Pi = list(res.np_data["pop_adi"][:, i, 0])
            plt.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
            plt.legend()
    
    
    plt.subplot(1, 2, 2)
    plt.title('Diabatic population dynamics' )
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    
    if "pop_dia" in properties_to_save and t != None:
        for i in range(nstates):
            #Pi = data_conv.unpack1(res.data["pop_dia"], i, 0, 0)                        
            Pi = list(res.np_data["pop_dia"][:, i, 0])
            plt.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
            plt.legend()
            
    plt.savefig("%s/Fig1.png" % (prefix), dpi=300)
    plt.savefig("%s/Fig1.pdf" % (prefix), dpi=300)
            
                
    #============= Energies =====================
    plt.figure(2, figsize=(36, 12)) # dpi=300, frameon=False)           
    
    plt.subplot(1, 2, 1)
    plt.title('Energies' )
    plt.xlabel('t, a.u.')
    plt.ylabel('Energy, a.u.')
    if "Ekin_dia" in properties_to_save \
       and "Epot_dia" in properties_to_save \
       and "Etot_dia" in properties_to_save \
       and t != None:
        
        #Ekin_dia = res.data["Ekin_dia"]  
        #Epot_dia = res.data["Epot_dia"]  
        #Etot_dia = res.data["Etot_dia"]  
        Ekin_dia = list(res.np_data["Ekin_dia"][:])
        Epot_dia = list(res.np_data["Epot_dia"][:])
        Etot_dia = list(res.np_data["Etot_dia"][:])

        plt.plot(t, Etot_dia, label='$Etot_{dia}$', linewidth=10, color = colors["11"])   
        plt.plot(t, Ekin_dia, label='$Ekin_{dia}$', linewidth=10, color = colors["21"])   
        plt.plot(t, Epot_dia, label='$Epot_{dia}$', linewidth=10, color = colors["31"])   
        plt.legend()
    
    plt.subplot(1, 2, 2)
    plt.title('Energies' )
    plt.xlabel('t, a.u.')
    plt.ylabel('Energy, a.u.')
    
    if "Ekin_adi" in properties_to_save \
       and "Epot_adi" in properties_to_save \
       and "Etot_adi" in properties_to_save \
       and t != None:
        
        #Ekin_adi = res.data["Ekin_adi"]  
        #Epot_adi = res.data["Epot_adi"]  
        #Etot_adi = res.data["Etot_adi"]  
        Ekin_adi = list(res.np_data["Ekin_adi"][:])
        Epot_adi = list(res.np_data["Epot_adi"][:])
        Etot_adi = list(res.np_data["Etot_adi"][:])

        plt.plot(t, Etot_adi, label='$Etot_{adi}$', linewidth=10, color = colors["11"])   
        plt.plot(t, Ekin_adi, label='$Ekin_{adi}$', linewidth=10, color = colors["21"])   
        plt.plot(t, Epot_adi, label='$Epot_{adi}$', linewidth=10, color = colors["31"])   
        plt.legend()

    plt.savefig("%s/Fig2.png" % (prefix), dpi=300)
    plt.savefig("%s/Fig2.pdf" % (prefix), dpi=300)
    
  
    #============= Phase spaces & Norms  =====================
    plt.figure(3, figsize=(36, 12)) # dpi=300, frameon=False)           
                
    plt.subplot(1, 2, 1)
    plt.title('Phase space' )
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Momentum, a.u.')   
    
    if "q_dia" in properties_to_save and "p_dia" in properties_to_save:
        #ndof = res.data["q_dia"][0].num_of_rows
        ndof = res.np_data["q_dia"].shape[1]
    
        for idof in range(ndof):
            #qi = data_conv.unpack1(res.data["q_dia"], idof, 0, 0)
            #pi = data_conv.unpack1(res.data["p_dia"], idof, 0, 0)
            qi = list(res.np_data["q_dia"][:, idof, 0])
            pi = list(res.np_data["p_dia"][:, idof, 0])
    
            plt.plot(qi, pi, label='', linewidth=10, color = colors[clrs_index[i]])   
            plt.legend()
                    

    plt.subplot(1, 2, 2)
    plt.title('Norms' )
    plt.xlabel('Time, a.u.')
    plt.ylabel('Norm')
    
    if "norm_dia" in properties_to_save and "norm_adi" in properties_to_save and t != None:        
        #nrm_dia = res.data["norm_adi"]
        #nrm_adi = res.data["norm_dia"]
        nrm_dia = list(res.np_data["norm_dia"][:])
        nrm_adi = list(res.np_data["norm_adi"][:])
                
        plt.plot(t, nrm_dia, label='Diabatic', linewidth=10, color = colors["11"])   
        plt.plot(t, nrm_adi, label='Adiabatic', linewidth=10, color = colors["21"])   
        plt.legend()

    plt.savefig("%s/Fig3.png" % (prefix), dpi=300)
    plt.savefig("%s/Fig3.pdf" % (prefix), dpi=300)
    
        
    plt.show()
    plt.close()



def plot_hdf5(plot_params, ax=plt, use_default_ax=True):
    """
    This function is meant to plot the results stored in the hdf files generated by the exact dynamics runs

    Args:

        prefix ( string ): the name of the directory containing the input HDF5 file
            This directory will also be used to output the generated picture files [ default : "out"]
        filename ( string ): name of the HDF5 file to read [ default: "data.hdf"]
        output_level ( int ): the level of info contained in the HDF5 file [ default : 3]        
        which_adi_states ( list of ints ) : indices of the adiabatic states to print [ default: [0] ]
        which_dia_states ( list of ints ) : indices of the diabatic states to print [ default: [0] ]
        colors ( dictionary ): the definition of the colors to use
        clrs_index ( list of strings ) : defines the mapping of the colors on integers and vice versa 
    

    """

    if use_default_ax==True:    
        ax.rc('axes', titlesize=24)      # fontsize of the axes title
        ax.rc('axes', labelsize=20)      # fontsize of the x and y labels
        ax.rc('legend', fontsize=20)     # legend fontsize
        ax.rc('xtick', labelsize=16)    # fontsize of the tick labels
        ax.rc('ytick', labelsize=16)    # fontsize of the tick labels

        ax.rc('figure.subplot', left=0.2)
        ax.rc('figure.subplot', right=0.95)
        ax.rc('figure.subplot', bottom=0.13)
        ax.rc('figure.subplot', top=0.88)

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
    default_params = {  "prefix":"out", "filename":"data.hdf", "hdf5_output_level":3,
                        "which_dofs":[0], "which_adi_states":[0], "which_dia_states":[0],
                        "colors":colors, "clrs_index":clrs_index,
                        "properties_to_save": 
                          [ "timestep", "time", "Ekin_dia", "Ekin_adi", "Epot_dia", 
                            "Epot_adi", "Etot_dia", "Etot_adi", "norm_dia", "norm_adi",
                            "pop_dia", "pop_adi", "q_dia", "q_adi", "q2_dia", "q2_adi", 
                            "p_dia", "p_adi", "p2_dia", "p2_adi",
                            "denmat_dia", "denmat_adi", 
                            "PSI_dia", "PSI_adi", "reciPSI_dia", "reciPSI_adi" ] 
                     }
    comn.check_input(plot_params, default_params, critical_params)
        

            
    filename = plot_params["filename"]
    prefix = plot_params["prefix"]
    hdf5_output_level = plot_params["hdf5_output_level"]    
    which_dofs = plot_params["which_dofs"]
    which_adi_states = plot_params["which_adi_states"]
    which_dia_states = plot_params["which_dia_states"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    properties_to_save = plot_params["properties_to_save"]

    out_prefix = prefix
    

    with h5py.File(F"{prefix}/{filename}", 'r') as f:
    
        t = None
        if "time" in properties_to_save:                        
            t = list(f["time/data"][:])
    
        #=============== Populations ======================
    
        ax.figure(1, figsize=(36, 12)) # dpi=300, frameon=False)    
        ax.subplot(1, 2, 1)
        ax.title('Adiabatic population dynamics' )
        ax.xlabel('Time, a.u.')
        ax.ylabel('Population')
    
        if "pop_adi" in properties_to_save and t != None:
            nstates = f["pop_adi/data"].shape[1] #.attrs['dim'][1]                        
            for i in range(nstates):        
                if i in which_adi_states:
                    Pi = list(f["pop_adi/data"][:, i, 0])                
                    ax.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
                    ax.legend()
    
    
        ax.subplot(1, 2, 2)
        ax.title('Diabatic population dynamics' )
        ax.xlabel('Time, a.u.')
        ax.ylabel('Population')
    
        if "pop_dia" in properties_to_save and t != None:
            nstates = f["pop_dia/data"].shape[1] #.attrs['dim'][1]                        
            for i in range(nstates):
                if i in which_dia_states:
                    Pi = list(f["pop_dia/data"][:, i, 0])                            
                    ax.plot(t, Pi, label='$P_%i$' % (i), linewidth=10, color = colors[clrs_index[i]])   
                    ax.legend()
            
        ax.savefig("%s/Fig1.png" % (prefix), dpi=300)
        ax.savefig("%s/Fig1.pdf" % (prefix), dpi=300)
            
                
        #============= Energies =====================
        ax.figure(2, figsize=(36, 12)) # dpi=300, frameon=False)           
    
        ax.subplot(1, 2, 1)
        ax.title('Energies' )
        ax.xlabel('t, a.u.')
        ax.ylabel('Energy, a.u.')
        if "Ekin_dia" in properties_to_save \
           and "Epot_dia" in properties_to_save \
           and "Etot_dia" in properties_to_save \
           and t != None:
        
            Ekin_dia =  list(f["Ekin_dia/data"][:])                             
            Epot_dia =  list(f["Epot_dia/data"][:])                             
            Etot_dia =  list(f["Etot_dia/data"][:])                             
        
            ax.plot(t, Etot_dia, label='$Etot_{dia}$', linewidth=10, color = colors["11"])   
            ax.plot(t, Ekin_dia, label='$Ekin_{dia}$', linewidth=10, color = colors["21"])   
            ax.plot(t, Epot_dia, label='$Epot_{dia}$', linewidth=10, color = colors["31"])   
            ax.legend()
    
    
        ax.subplot(1, 2, 2)
        ax.title('Energies' )
        ax.xlabel('t, a.u.')
        ax.ylabel('Energy, a.u.')
    
        if "Ekin_adi" in properties_to_save \
           and "Epot_adi" in properties_to_save \
           and "Etot_adi" in properties_to_save \
           and t != None:
        
            Ekin_adi =  list(f["Ekin_adi/data"][:])                             
            Epot_adi =  list(f["Epot_adi/data"][:])                             
            Etot_adi =  list(f["Etot_adi/data"][:])                             
        
            ax.plot(t, Etot_adi, label='$Etot_{adi}$', linewidth=10, color = colors["11"])   
            ax.plot(t, Ekin_adi, label='$Ekin_{adi}$', linewidth=10, color = colors["21"])   
            ax.plot(t, Epot_adi, label='$Epot_{adi}$', linewidth=10, color = colors["31"])   
            ax.legend()

        ax.savefig("%s/Fig2.png" % (prefix), dpi=300)
        ax.savefig("%s/Fig2.pdf" % (prefix), dpi=300)
    
  
        #============= Phase spaces & Norms  =====================
        ax.figure(3, figsize=(36, 12)) # dpi=300, frameon=False)           
                
        ax.subplot(1, 2, 1)
        ax.title('Phase space' )
        ax.xlabel('Coordinate, a.u.')
        ax.ylabel('Momentum, a.u.')   
    
        if "q_dia" in properties_to_save and "p_dia" in properties_to_save:                                
            ndof = f["q_dia/data"].shape[1] #.attrs['dim'][1]                        
    
            for idof in range(ndof):
                if idof in which_dofs:
                    qi =  list(f["q_dia/data"][:, idof, 0])                             
                    pi =  list(f["p_dia/data"][:, idof, 0])                                         
    
                    ax.plot(qi, pi, label='', linewidth=10, color = colors[clrs_index[i]])   
                    ax.legend()
                    
                
        ax.subplot(1, 2, 2)
        ax.title('Norms' )
        ax.xlabel('Time, a.u.')
        ax.ylabel('Norm')
    
        if "norm_dia" in properties_to_save and "norm_adi" in properties_to_save and t != None:
        
            nrm_dia = list(f["norm_dia/data"][:])                             
            nrm_adi = list(f["norm_adi/data"][:])                                     
                
            ax.plot(t, nrm_dia, label='Diabatic', linewidth=10, color = colors["11"])   
            ax.plot(t, nrm_adi, label='Adiabatic', linewidth=10, color = colors["21"])   
            ax.legend()

        ax.savefig("%s/Fig3.png" % (prefix), dpi=300)
        ax.savefig("%s/Fig3.pdf" % (prefix), dpi=300)
            
        ax.show()
        ax.close()
        
