#*********************************************************************************                     
#* Copyright (C) 2019-2020 Brendan A. Smith, Alexey V. Akimov
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
   :synopsis: This module implements the functions to read the HDF5 files and plot various properties
       List of functions:
           * plot_dyn(plot_params)
           * plot_dyn_old(res)

.. moduleauthor:: Alexey V. Akimov
  
"""

__author__ = "Brendan A. Smith, Alexey V. Akimov"
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
import matplotlib.pyplot as plt 
#from matplotlib.mlab import griddata

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import libra_py.data_visualize as davis
import libra_py.units as units

import util.libutil as comn



def common_defaults(plot_params_):
    """
    This function returns a dictionary of the plotting parameters
    set up as either the input values provided by the user (from `plot_params_`), if 
    given, or the intrinsic defaults values, if not defined by the user
    
    Args:
        plot_params_ ( dict ): plotting parameters dictionary, may contain the following keys:
        
            * **prefix** ( string ): directory containing the data files, it also serves 
                as the directory where to save the figures, if requested [ default: "out" ]
            
            * **filename** ( string ): the name of the HDF data file containing the information to 
                be plotted [ default: "mem_data.hdf"]
                
            * **output_level** ( int ): the level of info contained in the HDF5 file [ default : 3 ]
            
            * **which_trajectories** ( list of ints ): indices of the trajectories to plot [ default: [0] ]
            
            * **which_dofs** ( list of ints ): indices of the DOFs to plot [ default: [0] ]
            
            * **which_adi_states** ( list of ints ) : indices of the adiabatic states to plot [ default: [0] ]
            
            * **which_dia_states** ( list of ints ) : indices of the diabatic states to print [ default: [0] ]
            
            * **which_energies** (list of strings ): determine which types of energies to plot 
                Can contain any of the following: 
                
                - "kinetic": kinetic energy
                - "potential": potential energy (of the active state)
                - "total": kinetic + potential - shall be conserved in NVE
                - "extended": total + bath (thermostat) - shall be conserved in NVT
                [ default: ["kinetic", "potential", "total"] ]
            
            * **colors** ( dictionary ): the definition of the colors to use
            
            * **clrs_index** ( list of strings ) : defines the mapping of the colors on integers 
                and vice versa 

            * **save_figures** ( int ): whether to save the figures as png files in the `prefix` directory  
                - 0: no, but the figures will still be shown in Jupyter notebook and can be saved manually
                - 1: yes [ default: 1]
                                            
            * **dpi** ( int ): the DPI level to control the quality of the figures produces [ default: 300 ]
            
            * **figsize** (2-double tuple): define the size of the figures in inches [ default: (6.42, 2.41) ]
            
            * **axes_fontsize** ( tuple of 2 ints ): font sizes for x and y tics [ default: (8, 8) ]
            
            * **axes_label_fontsize** ( tuple of 2 ints ): font sizes for x and y axes labels [ default: (8, 8) ]
            
            * **legend_fontsize** ( int ): legend font size [ default: 8 ]
            
            * **title_fontsize** ( int ): title font size [ default: 8 ]
            
            * **linewidth** ( int ): linewidth [ default: 2 ]
            
            * **frameon** ( boolean ): whether to use Frame on each panel [ default: True ]
            
            * **xlim** ( tuple of 2 doubles or None): 
                None: x boundaries of the plot are determined automatically [ default: None ]
                tuple of 2 doubles: use them to determine the range of the x axis
                
            * **ylim** ( tuple of 2 doubles or None): 
                None: y boundaries of the plot are determined automatically [ default: None ]
                tuple of 2 doubles: use them to determine the range of the y axis 
    
    """
    
    
    plot_params = dict(plot_params_ )
    
    
    critical_params = [  ] 
    default_params = {  "prefix":"out", "filename":"mem_data.hdf", "output_level":3,
                        "which_trajectories":[0], "which_dofs":[0],
                        "which_adi_states":[0], "which_dia_states":[0],
                        "which_energies":["kinetic", "potential", "total"],
                        "colors":davis.colors, "clrs_index":davis.clrs_index,
                        "save_figures":1, "do_show":1, "dpi":300, "figsize":(6.42, 2.41),
                        "axes_fontsize":(8, 8), "axes_label_fontsize":(8,8),
                        "legend_fontsize":8, "title_fontsize":8,
                        "linewidth":2,
                        "frameon":True,
                        "xlim":None, "ylim":None,
                        "what_to_plot":["energies", "energy_fluctuations", 
                                        "traj_resolved_adiabatic_ham", "traj_resolved_diabatic_ham",
                                        "coordinates", "phase_space",
                                        "sh_pop", "se_pop_adi", "se_pop_dia",
                                        "sh_pop_raw", "se_pop_adi_raw", "se_pop_dia_raw",
                                        "time_overlaps", "projectors", "basis_transform"
                                       ],
                        "no_label":0                      
                     }
    comn.check_input(plot_params, default_params, critical_params)

    
    return plot_params


def add_energies(plt, hdf_file, plot_params_, property_type):
    """
    Adds the plotting of the ensemble-averaged kinetic, potential, and total energies
    and their fluctuations
    This function does not plot it though
    
    property_type ( string): what kind of data to plot, can be one of ["energies", "energy_fluctuations"]    
    
    Call it after adding a sub-plot
    """

    possible_options = ["energies", "energy_fluctuations"]
    if property_type not in possible_options:
        print(F"Error in add_energies - the property_type argument {property_type} is invalid\n")
        print(F"Must be one of the following options: {possible_options}\nExiting")
        sys.exit(0)

    
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]   
    which_energies = plot_params["which_energies"]

            
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
    titles = { "energies": "Ensemble Averaged Energies",
               "energy_fluctuations": "Energy Fluctuations"               
             }


    plt.title(titles[property_type], fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Energy, a.u.', fontsize=axes_label_fontsize[1])        
        

    res = 0

    if property_type == "energies":        
        if "potential" in which_energies and "Epot_ave/data" in hdf_file.keys() and "time/data" in hdf_file.keys():
            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["Epot_ave/data"][:],
                     label="Potential energy", linewidth=Lw, color = colors[ clrs_index[0] ])
            res = 1
        if "kinetic" in which_energies and "Ekin_ave/data" in hdf_file.keys() and "time/data" in hdf_file.keys():
            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["Ekin_ave/data"][:],
                     label="Kinetic energy", linewidth=Lw, color = colors[ clrs_index[1] ])
            res = 1
        if "total" in which_energies and "Etot_ave/data" in hdf_file.keys() and "time/data" in hdf_file.keys():
            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["Etot_ave/data"][:], 
                     label="Total energy", linewidth=Lw, color = colors[ clrs_index[2] ])
            res = 1
        if "extended" in which_energies and "E_NHC/data" in hdf_file.keys() and "time/data" in hdf_file.keys():
            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["E_NHC/data"][:],
                     label="Extended energy", linewidth=Lw, color = colors[ clrs_index[3] ])
            res = 1
                        
    elif property_type == "energy_fluctuations":
        if "potential" in which_energies and "dEpot_ave/data" in hdf_file.keys() and "time/data" in hdf_file.keys():
            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["dEpot_ave/data"][:],
                     label="Potential energy", linewidth=Lw, color = colors[ clrs_index[0] ])
            res = 1
        if "kinetic" in which_energies and "dEkin_ave/data" in hdf_file.keys() and "time/data" in hdf_file.keys():
            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["dEkin_ave/data"][:],
                     label="Kinetic energy", linewidth=Lw, color = colors[ clrs_index[1] ])
            res = 1
        if "total" in which_energies and "dEtot_ave/data" in hdf_file.keys() and "time/data" in hdf_file.keys():
            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["dEtot_ave/data"][:],
                     label="Total energy", linewidth=Lw, color = colors[ clrs_index[2] ])
            res = 1
        # This one is not yet present
        #if "extended" in which_energies:
        #    plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["dE_NHC/data"][:],
        #             label="Extended energy", linewidth=Lw, color = colors[ clrs_index[3] ])                    
    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()

    return res
                

def add_trajectory_resolved_ham_property(plt, hdf_file, plot_params_, ham_property_type):
    """
    Adds the plotting of the trajectory-resolved Hamiltonian-related properties
    This function does not plot it though
    
    Args: 
    
        ham_property_type ( string ): selects what kind of property to plot, can be one of the following:
        
            - hvib_adi: diagonal elements of the adiabatic Hamiltonian
            - hvib_dia: diagonal elements of the diabatic Hamiltonian        
    
    Call it after adding a sub-plot
    """
    
    possible_options = ["hvib_adi", "hvib_dia"]
    if ham_property_type not in possible_options:
        print(F"Error in add_trajectory_resolved_ham_property - the ham_property_type argument\
              {ham_property_type} is invalid\n")
        print(F"Must be one of the following options: {possible_options}\nExiting")
        sys.exit(0)
        
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]   
            
    which_trajectories = plot_params["which_trajectories"]
    which_states, nstates, ntraj = [], 0, 0
    no_label = plot_params["no_label"]


    if ham_property_type == "hvib_adi":
        which_states = plot_params["which_adi_states"]
        ntraj = hdf_file["hvib_adi/data"].shape[1] 
        nstates = hdf_file["hvib_adi/data"].shape[2]
    
    elif ham_property_type == "hvib_dia":
        which_states = plot_params["which_dia_states"]
        ntraj = hdf_file["hvib_dia/data"].shape[1] 
        nstates = hdf_file["hvib_dia/data"].shape[2]
    
        
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
    labels = { "hvib_adi": "Adiabatic state energies",
               "hvib_dia": "Diabatic state energies"        
             }
            
    plt.title(labels[ham_property_type], fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Energy, a.u.', fontsize=axes_label_fontsize[1])   
    


    res = 0 
    indx = -1
    for istate in range(nstates):
        if istate in which_states:
            indx = indx + 1
            for tr in range(ntraj):
                if tr in which_trajectories:                    
                    if "time/data" in hdf_file.keys() and F"{ham_property_type}/data" in hdf_file.keys():
                        if no_label:
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, 
                                     hdf_file[F"{ham_property_type}/data"][:, tr, istate, istate],
                                     label="",
                                     linewidth=Lw, color = colors[ clrs_index[indx] ])
                            res = 1

                        else:
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, 
                                     hdf_file[F"{ham_property_type}/data"][:, tr, istate, istate],
                                     label=F"traj={tr}, state={istate}",
                                     linewidth=Lw, color = colors[ clrs_index[indx] ])
                            res = 1

    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()

    return res
    

def add_cooordinates_vs_t(plt, hdf_file, plot_params_):
    """
    Adds the plotting of the trajectory and dof-resolved coordinates vs. time
    This function does not plot it though
    
    Call it after adding a sub-plot
    """
    
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]   
    
    which_trajectories = plot_params["which_trajectories"]
    which_dofs   = plot_params["which_dofs"]

    no_label = plot_params["no_label"]
    
        
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
    
    plt.title('Time-dependent Coordinates, a.u.', fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Coordinate, a.u.', fontsize=axes_label_fontsize[1])   



    res = 0
    ntraj, ndof = 0, 0    

    if "q/data" in hdf_file.keys():
        ntraj = hdf_file["q/data"].shape[1]
        ndofs = hdf_file["q/data"].shape[2]
        res = 1

    if res==1 and "time/data" in hdf_file.keys():    
        for tr in range(ntraj):        
            if tr in which_trajectories:
                for dof in range(ndofs):
                    if dof in which_dofs:
                        if no_label:
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["q/data"][:, tr, dof],
                                     label="", linewidth=Lw, 
                                     color = colors[ clrs_index[dof] ]) 
                        else: 
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["q/data"][:, tr, dof],
                                     label=F"traj={tr} dof={dof}", linewidth=Lw, 
                                     color = colors[ clrs_index[dof] ]) 
    else:
        res = 0
    
    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()

    return res

        

def add_momenta_vs_t(plt, hdf_file, plot_params_):
    """
    Adds the plotting of the trajectory and dof-resolved momenta vs. time
    This function does not plot it though
    
    Call it after adding a sub-plot
    """
    
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]   
    
    which_trajectories = plot_params["which_trajectories"]
    which_dofs   = plot_params["which_dofs"]

    no_label = plot_params["no_label"]

            
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
    
    plt.title('Time-dependent Momenta, a.u.', fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Momentum, a.u.', fontsize=axes_label_fontsize[1])   


    res = 0
    ntraj, ndof = 0, 0    

    if "p/data" in hdf_file.keys():    
        ntraj = hdf_file["p/data"].shape[1]
        ndofs = hdf_file["p/data"].shape[2]
        res = 1

    if res==1 and "time/data" in hdf_file.keys():        
        for tr in range(ntraj):        
            if tr in which_trajectories:
                for dof in range(ndofs):
                    if dof in which_dofs:
                        if no_label:
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["p/data"][:, tr, dof],
                                     label="", linewidth=Lw, 
                                     color = colors[ clrs_index[dof] ]) 
                        else:
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["p/data"][:, tr, dof],
                                     label=F"traj={tr} dof={dof}", linewidth=Lw, 
                                     color = colors[ clrs_index[dof] ]) 
    else:
        res =  0

    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()

    return res        


def add_forces_vs_t(plt, hdf_file, plot_params_):
    """
    Adds the plotting of the trajectory and dof-resolved momenta vs. time
    This function does not plot it though

    Call it after adding a sub-plot
    """

    plot_params = common_defaults(plot_params_)

    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]

    which_trajectories = plot_params["which_trajectories"]
    which_dofs   = plot_params["which_dofs"]

    no_label = plot_params["no_label"]


    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
    plt.title('Time-dependent Forces, a.u.', fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Force, a.u.', fontsize=axes_label_fontsize[1])


    res = 0
    ntraj, ndof = 0, 0

    if "f/data" in hdf_file.keys():
        ntraj = hdf_file["f/data"].shape[1]
        ndofs = hdf_file["f/data"].shape[2]
        res = 1

    if res==1 and "time/data" in hdf_file.keys():
        for tr in range(ntraj):
            if tr in which_trajectories:
                for dof in range(ndofs):
                    if dof in which_dofs:
                        if no_label:
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["f/data"][:, tr, dof],
                                     label="", linewidth=Lw,
                                     color = colors[ clrs_index[dof] ])
                        else:
                            plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file["f/data"][:, tr, dof],
                                     label=F"traj={tr} dof={dof}", linewidth=Lw,
                                     color = colors[ clrs_index[dof] ])
    else:
        res =  0

    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()

    return res



def add_phase_space(plt, hdf_file, plot_params_):
    """
    Adds the plotting of the trajectory and dof-resolved momenta vs. coordinates (phase space)
    This function does not plot it though
    
    Call it after adding a sub-plot
    """
    
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]   
    
    which_trajectories = plot_params["which_trajectories"]
    which_dofs   = plot_params["which_dofs"]

    no_label = plot_params["no_label"]

            
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        

    plt.title('Phase Space Portrait', fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Coordinate, a.u.', fontsize=axes_label_fontsize[0])
    plt.ylabel('Momentum, a.u.', fontsize=axes_label_fontsize[1])   

    res = 0
    ntraj, ndof = 0, 0    
    if "q/data" in hdf_file.keys() and "p/data" in hdf_file.keys() :    
        ntraj = hdf_file["p/data"].shape[1]
        ndofs = hdf_file["p/data"].shape[2]
        res = 1

    if res==1:
        for tr in range(ntraj):        
            if tr in which_trajectories:
                for dof in range(ndofs):
                    if dof in which_dofs:
                        if no_label:
                            plt.plot(hdf_file["q/data"][:, tr, dof], hdf_file["p/data"][:, tr, dof],
                                     label="", linewidth=Lw, 
                                     color = colors[ clrs_index[dof] ]) 
                        else:
                            plt.plot(hdf_file["q/data"][:, tr, dof], hdf_file["p/data"][:, tr, dof],
                                     label=F"traj={tr} dof={dof}", linewidth=Lw, 
                                     color = colors[ clrs_index[dof] ]) 
        
    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()
                
    return res    
            
        
def add_populations(plt, hdf_file, plot_params_, pop_type ):
    """
    Adds the plotting of the ensemble-averaged SE or SH populations vs. time
    This function does not plot it though
    
    
    pop_type (string): selector of the type of the population to plot
    
        - "D_dia_raw"
        - "D_adi_raw"
        - "SH_pop_raw"
        - "D_dia"
        - "D_adi"
        - "SH_pop"
        
    
    Call it after adding a sub-plot
    """
    
    possible_options = ["D_dia_raw", "D_adi_raw", "SH_pop_raw", "D_dia", "D_adi", "SH_pop", "se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia"]
    if pop_type not in possible_options:
        print(F"Error in add_populations - the pop_type argument {pop_type} is invalid\n")
        print(F"Must be one of the following options: {possible_options}\nExiting")
        sys.exit(0)            
    
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"] * 10  
    
    
    nstates, which_states = 0, []

    if pop_type in ["D_dia"]: #, "D_dia_raw"]:  # diabatic SE properties 
        if "D_dia/data" in hdf_file.keys():          
            nstates = hdf_file["D_dia/data"].shape[1] 
            which_states = plot_params["which_dia_states"]

    elif pop_type in ["se_pop_dia"]:  # diabatic SE populations
        if "se_pop_dia/data" in hdf_file.keys():          
            nstates = hdf_file["se_pop_dia/data"].shape[1] 
            which_states = plot_params["which_dia_states"]

        
    elif pop_type in ["D_adi"]: #, "D_adi_raw"]: # adiabatic SE properties
        if "D_adi/data" in hdf_file.keys():
            nstates = hdf_file["D_adi/data"].shape[1]
            which_states = plot_params["which_adi_states"]

    elif pop_type in ["se_pop_adi"]:  # adiabatic SE populations
        if "se_pop_adi/data" in hdf_file.keys():          
            nstates = hdf_file["se_pop_adi/data"].shape[1] 
            which_states = plot_params["which_adi_states"]

        
    elif pop_type in ["SH_pop"]: #, "SH_pop_raw"]: # adiabatic SH properties
        if "SH_pop/data" in hdf_file.keys():
            nstates = hdf_file["SH_pop/data"].shape[1]
            which_states = plot_params["which_adi_states"]

    elif pop_type in ["sh_pop_dia"]:  # diabatic SH populations
        if "sh_pop_dia/data" in hdf_file.keys():          
            nstates = hdf_file["sh_pop_dia/data"].shape[1] 
            which_states = plot_params["which_dia_states"]

    elif pop_type in ["sh_pop_adi"]:  # adiabatic SH populations
        if "sh_pop_adi/data" in hdf_file.keys():          
            nstates = hdf_file["sh_pop_adi/data"].shape[1] 
            which_states = plot_params["which_adi_states"]

        
    titles = { "D_dia": "Diabatic SE populations",
               "D_dia_raw": "Diabatic SE populations (raw)",
               "se_pop_dia": "Diabatic SE populations",
               "D_adi": "Adiabatic SE populations",
               "D_adi_raw": "Adiabatic SE populations (raw)",
               "se_pop_adi": "Adiabatic SE populations",
               "SH_pop": "Adiabatic SH populations",
               "SH_pop_raw": "Adiabatic SH populations (raw)",
               "sh_pop_dia": "Diabatic SH populations",
               "sh_pop_adi": "Adiabatic SH populations"
             }
    
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
            
    plt.title(titles[pop_type], fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Population', fontsize=axes_label_fontsize[1])   
            
#    if pop_type in ["D_dia_raw", "D_adi_raw", "D_dia", "D_adi"]:

    res = 0
    if pop_type in ["D_dia", "D_adi"]:
        if F"{pop_type}/data" in hdf_file.keys():
            res = 1
            indx = -1
            for istate in range(nstates):
                if istate in which_states:
                    indx = indx + 1            
                    plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file[F"{pop_type}/data"][:, istate, istate], 
                             label=F"state {istate}", linewidth=Lw, color = colors[clrs_index[indx] ])
                
    elif pop_type in ["SH_pop"]:
        if F"{pop_type}/data" in hdf_file.keys():
            res = 1
            indx = -1
            for istate in range(nstates):
                if istate in which_states:
                    indx = indx + 1            
                    plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file[F"{pop_type}/data"][:, istate, 0], 
                             label=F"state {istate}", linewidth=Lw, color = colors[clrs_index[indx] ]) 


    elif pop_type in ["se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia"]:
        if F"{pop_type}/data" in hdf_file.keys():
            res = 1
            indx = -1
            for istate in range(nstates):
                if istate in which_states:
                    indx = indx + 1            
                    plt.plot(hdf_file["time/data"][:]/units.fs2au, hdf_file[F"{pop_type}/data"][:, istate], 
                             label=F"state {istate}", linewidth=Lw, color = colors[clrs_index[indx] ]) 

        
    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()            
        
    return res



def add_time_overlaps_projectors(plt, hdf_file, plot_params_, prop_type):
    """
    Adds the plotting of the time-overlaps or projectors vs. time    
    This function does not plot it though
        
    prop_type (string): selector of the type of the property to plot
    
        - "St" : <psi_i(t-dt)|psi_i(t)>, for selected adiabatic states,
            this property indicates the state-reorderings
        - "projector": the diagonal elements of the cumulative transformation matrix
            that converts the "raw" adiabatic states to the dynamically-consistent ones                            
    
    Call it after adding a sub-plot
    """
    
    possible_options = ["St", "projector"]
    if prop_type not in possible_options:
        print(F"Error in add_time_overlaps_projectors - the prop_type argument {prop_type} is invalid\n")
        print(F"Must be one of the following options: {possible_options}\nExiting")
        sys.exit(0)            
    
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]   
    
    which_states = plot_params["which_adi_states"]
    which_trajectories = plot_params["which_trajectories"]
    
    ntraj = hdf_file[F"{prop_type}/data"].shape[1] 
    nstates = hdf_file[F"{prop_type}/data"].shape[2] 
                    
    titles = { "St": "Time-overlaps of adiabatic states",
               "projector": "Projectors of raw to dyn-consistent states",               
             }
    
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
            
    plt.title(titles[prop_type], fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Magnitude', fontsize=axes_label_fontsize[1])   

    
    res = 0
    indx = -1
    for istate in range(nstates):
        if istate in which_states:
            indx = indx + 1
            for tr in range(ntraj):
                if tr in which_trajectories:
                    if "time/data" in hdf_file.keys() and F"{prop_type}/data" in hdf_file.keys():
                        plt.plot(hdf_file["time/data"][:]/units.fs2au, 
                                 hdf_file[F"{prop_type}/data"][:, tr, istate, istate], 
                                 label=F"state {istate}", linewidth=Lw, color = colors[clrs_index[indx] ])
                        res = 1
                    
    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()   
    
    return res
    

def add_basis_transform(plt, hdf_file, plot_params_):
    """
    Adds the plotting of the adiabatic-to-diabatic transforms vs. time    
    This function does not plot it though    
    
    Call it after adding a sub-plot
    """
    
    plot_params = common_defaults(plot_params_)
            
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]   
    
    which_trajectories = plot_params["which_trajectories"]
    which_dia_states = plot_params["which_dia_states"]
    which_adi_states = plot_params["which_adi_states"]
        
    ntraj = hdf_file["basis_transform/data"].shape[1]
    ndia  = hdf_file["basis_transform/data"].shape[2]
    nadi  = hdf_file["basis_transform/data"].shape[3]
                    
        
    if xlim!=None:
        plt.xlim( xlim[0], xlim[1])
    if ylim!=None:
        plt.ylim( ylim[0], ylim[1])
        
        
    plt.title("Adiabatic to diabatic transformation", fontsize=title_fontsize)
    plt.xticks(fontsize=axes_fontsize[0])
    plt.yticks(fontsize=axes_fontsize[1])                            
    plt.xlabel('Time, fs', fontsize=axes_label_fontsize[0])
    plt.ylabel('Magnitude', fontsize=axes_label_fontsize[1])   
            
    res = 0 
    indx = -1
    for istate in range(nadi):
        if istate in which_adi_states:
            for istate2 in range(ndia):
                if istate2 in which_dia_states:
                            
                    indx = indx + 1
                    
                    lbl='$< \psi^{dia}_{%i} | \psi^{adi}_{%i} >$' % (istate2, istate)
                    
                    for tr in range(ntraj):
                        if tr in which_trajectories:                                     
                            if "time/data" in hdf_file.keys() and "basis_transform/data" in hdf_file.keys():
                                plt.plot(hdf_file["time/data"][:]/units.fs2au,
                                         hdf_file["basis_transform/data"][:, tr, istate2, istate], 
                                         label=lbl, linewidth=Lw, 
                                         color = colors[ clrs_index[indx] ])     
                                res = 1
                            
    plt.legend(fontsize=legend_fontsize)
    plt.tight_layout()   

    return res        
        
    
def plot_dynamics(plot_params_):
    """
    This function plots a 2 x 1 figure containing:
    supblot 1: average kinetic, potential, and total energies for the ensemble of trajectories
    subplot 2: the adiabatic potential energies for selected trajectories
        
    """
    
    plot_params = common_defaults(plot_params_)
        
    
    filename           = plot_params["filename"]
    prefix             = plot_params["prefix"]
    output_level       = plot_params["output_level"]
    which_dofs         = plot_params["which_dofs"]
    which_trajectories = plot_params["which_trajectories"]
    which_adi_states   = plot_params["which_adi_states"]
    which_dia_states   = plot_params["which_dia_states"]
    what_to_plot       = plot_params["what_to_plot"]  
    out_prefix = prefix
    
    axes_fontsize = plot_params["axes_fontsize"]
    axes_label_fontsize = plot_params["axes_label_fontsize"]
    legend_fontsize = plot_params["legend_fontsize"]
    title_fontsize = plot_params["title_fontsize"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"]
    Lw = plot_params["linewidth"]
    do_show = plot_params["do_show"]    

    with h5py.File(F"{prefix}/{filename}", 'r') as f:
        
        #===== Energy Components vs. Time =========
        if "energies" in what_to_plot:
            plt.figure(num=1, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                   edgecolor='black', frameon=plot_params["frameon"])        
            plt.subplot(1,1,1)                
            res = add_energies(plt, f, plot_params, "energies")        
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/average_energies-vs-t.png", dpi=plot_params["dpi"])


        #===== Energy Components Fluctuations vs. Time =========
        if "energy_fluctuations" in what_to_plot:
            plt.figure(num=2, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                   edgecolor='black', frameon=plot_params["frameon"])        
            plt.subplot(1,1,1)
            res = add_energies(plt, f, plot_params, "energy_fluctuations")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/average_energy_fluctuations-vs-t.png", dpi=plot_params["dpi"])


        #===== Coordinates vs Time  =========            
        if "coordinates" in what_to_plot:
            plt.figure(num=3, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_cooordinates_vs_t(plt, f, plot_params)
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/q-vs-t.png", dpi=plot_params["dpi"])


        #===== Momenta vs Time  =========            
        if "momenta" in what_to_plot:
            plt.figure(num=4, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_momenta_vs_t(plt, f, plot_params)
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/p-vs-t.png", dpi=plot_params["dpi"])

        #===== Forces vs Time  =========
        if "forces" in what_to_plot:
            plt.figure(num=5, figsize=plot_params["figsize"], dpi=plot_params["dpi"],
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_forces_vs_t(plt, f, plot_params)
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/f-vs-t.png", dpi=plot_params["dpi"])

        #===== Phase space  =========            
        if "phase_space" in what_to_plot:
            plt.figure(num=6, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_phase_space(plt, f, plot_params)
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/phase_space.png", dpi=plot_params["dpi"])


        #====== Populations of all kinds ================
        if "se_pop_adi" in what_to_plot:
            plt.figure(num=7, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])        
            plt.subplot(1,1,1)            
            res = add_populations(plt, f, plot_params_, "se_pop_adi")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/se_pop_adi.png", dpi=plot_params["dpi"])

        if "se_pop_dia" in what_to_plot:
            plt.figure(num=8, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])        
            plt.subplot(1,1,1)            
            res = add_populations(plt, f, plot_params_, "se_pop_dia")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/se_pop_dia.png", dpi=plot_params["dpi"])

        if "sh_pop_adi" in what_to_plot:
            plt.figure(num=9, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])        
            plt.subplot(1,1,1)            
            res = add_populations(plt, f, plot_params_, "sh_pop_adi")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/sh_pop_adi.png", dpi=plot_params["dpi"])

        if "sh_pop_dia" in what_to_plot:
            plt.figure(num=10, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])        
            plt.subplot(1,1,1)            
            res = add_populations(plt, f, plot_params_, "sh_pop_dia")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/sh_pop_dia.png", dpi=plot_params["dpi"])
                                       
                       
        #===== Trajectory-resolved adiabatic energies =========
        if "traj_resolved_adiabatic_ham" in what_to_plot:
            plt.figure(num=11, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_trajectory_resolved_ham_property(plt, f, plot_params, "hvib_adi")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/hvib_adi.png", dpi=plot_params["dpi"])

        if "traj_resolved_adiabatic_ham" in what_to_plot:
            plt.figure(num=12, figsize=plot_params["figsize"], dpi=plot_params["dpi"],
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_trajectory_resolved_ham_property(plt, f, plot_params, "hvib_dia")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/hvib_dia.png", dpi=plot_params["dpi"])


        #===== Time-overlaps and projectors =========
        if "time_overlaps" in what_to_plot:
            plt.figure(num=13, figsize=plot_params["figsize"], dpi=plot_params["dpi"],
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_time_overlaps_projectors(plt, f, plot_params, "St")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/time_overlaps.png", dpi=plot_params["dpi"])


        if "projector" in what_to_plot:
            plt.figure(num=14, figsize=plot_params["figsize"], dpi=plot_params["dpi"],
                       edgecolor='black', frameon=plot_params["frameon"])
            plt.subplot(1,1,1)
            res = add_time_overlaps_projectors(plt, f, plot_params, "projector")
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/projectors.png", dpi=plot_params["dpi"])


        #===== Basis transforms =========
        if "basis_transform" in what_to_plot:
            plt.figure(num=15, figsize=plot_params["figsize"], dpi=plot_params["dpi"], 
                       edgecolor='black', frameon=plot_params["frameon"])          
            plt.subplot(1,1,1)
            res = add_basis_transform(plt, f, plot_params)                   
            if plot_params["save_figures"]==1 and res==1:
                plt.savefig(F"{out_prefix}/basist_transform-vs-t.png", dpi=plot_params["dpi"])            

        
        if do_show:         
            plt.show()



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
          
            #========= SH, adi SE and dia SE populations =========
            nadi = f["D_adi_raw/data"].shape[1] #attrs['dim'][1]
            ndia = f["D_dia_raw/data"].shape[1] #attrs['dim'][1]

            plt.figure(num=None, figsize=(24, 8), dpi=300, frameon=False)        

            #================ SH populations =============
            plt.subplot(1,3,1)            
            plt.title('SH populations (raw)')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Population, a.u.')        
                        
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    indx = indx + 1
                    plt.plot(f["time/data"][:], f["SH_pop_raw/data"][:, istate, 0], label=F"state {istate}", linewidth=2, color = colors[clrs_index[indx] ])                 

            #================ adi SE populations =============
            plt.subplot(1,3,2)            
            plt.title('adi SE populations (raw)')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Population, a.u.')        
            
            indx = -1
            for istate in range(nadi):
                if istate in which_adi_states:
                    indx = indx + 1
                    plt.plot(f["time/data"][:], f["D_adi_raw/data"][:, istate, istate], label=F"state {istate}", linewidth=2, color = colors[clrs_index[indx] ])                 
                    
                    
            #================ dia SE populations =============
            plt.subplot(1,3,3)            
            plt.title('dia SE populations (raw)')
            plt.xlabel('Time, a.u.')
            plt.ylabel('Population, a.u.')        
            
            indx = -1
            for istate in range(ndia):
                if istate in which_dia_states:
                    indx = indx + 1
                    plt.plot(f["time/data"][:], f["D_dia_raw/data"][:, istate, istate], label=F"state {istate}", linewidth=2, color = colors[clrs_index[indx] ])                 
                    
                    
            plt.savefig(F"{out_prefix}/t-pops_raw.png", dpi=300)
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



def hdf2xyz(labels, filename, snaps, trajectories, atoms, unit_conversion_factor=1.0):
    """
    This function creates a string containing an xyz-formatted trajectory (multiple geometries)
    from an HDF5 file produced by the `compute.run_dynamics` function

    Args:
        labels ( list of strings ): the labels of the atoms of the system, the labels are in a 
            specific order that corresponds to the atomic order in the system. The length of this 
            list should be `nat` - the number of atoms listed in the xyz header, even if not all of 
            them are listed in the input parameter `atoms`
        filename ( string ): the name of the HDF5 file that contains the geometry 
        snaps ( list of ints ): indices of the timesteps to be included in the xyz file being generated.
            This allows printing only simesteps of interest. 
        trajectories ( list of ints ): indices of the trajectories that we want to include in the xyz file. 
            This allows plotting many trajectories at once (in every frame).
        atoms ( list of ints ): indices of the atomic species to be printed out to the xyz files (e.g. to 
            eventually visualize). This allows for visualizing only a subset of atoms (e.g. a spatial region,
            a molecule, a group, etc.) of the main interest.
        unit_conversion_factor ( double ): the conversion factor to convert the data stored in the
            HDF5 file to the new units. Since `compute.run_dynamics` function stores all the information
            in the atomic units (Bohrs for the coordinates) and since the xyz files visualization works
            best in the Angstrom units, it is common situation to use the 
            `unit_conversion_factor = 1.0/units.Angst` conversion factor.
 
    Returns:
        string: the string representation of the xyz trajectory file that is made of the 
            provided input geomtries/trajectories/atomic labels 
    
    """
    
    natoms = len(atoms)  # the actual number of atoms to show 
    ntraj = len(trajectories)
    
    md_xyz = ""
    
    with h5py.File(filename, 'r') as f:
    
        for isnap in snaps:
            
            md_xyz = md_xyz + F"{natoms*ntraj}\nsnapshot {isnap}\n"
            
            for itraj in trajectories:
                for iatom in atoms:
                
                    x = f["q/data"][isnap, itraj, 3*iatom+0] * unit_conversion_factor
                    y = f["q/data"][isnap, itraj, 3*iatom+1] * unit_conversion_factor
                    z = f["q/data"][isnap, itraj, 3*iatom+2] * unit_conversion_factor
    
                    md_xyz = md_xyz + F"{labels[iatom] }  {x} {y} {z}\n"                    
    
    return md_xyz
