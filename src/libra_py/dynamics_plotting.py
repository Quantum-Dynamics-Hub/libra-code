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


def plot_pes_properties(comp_model, model_params, pes_params_, plot_params_):    
    """
        Args:
        
            comp_model ( PyObject ): the function that returns the class with Hamiltonian properties
            
            model_params ( dictionary ): parameters of the model Hamiltonian
            
            pes_params_ ( dictionary ): controls the way the calculations are done
            
                Can contain the following parameters
                
                * **pes_params["ndia"]** ( int ): dinemsionality of the diabatic Hamiltonian [ default: 2 ]
                * **pes_params["nadi"]** ( int ): dinemsionality of the adiabatic Hamiltonian [ default: 2 ]
                * **pes_params["ndof"]** ( int ): the number of nuclear DOF [ default: 1 ]
                * **pes_params["active_dof"]** ( int ): index of the DOF w.r.t. which the scan will be done [ default: 0 ]
                * **pes_params["coord_type"]** ( 0 or 1 ): 0 - actual coordinate, 1 - time coordinate [ default: 0 ]
                * **pes_params["reference_coord"]** ( MARIX(ndof, 1) ): the geometry of a many-DOF-system, used in case we 
                    only vary one DOF, but want to keep all other DOFs fixed
                * **pes_params["coord_mapping"]** ( PyObject ): the function that maps the 1D "scan coordinate" to a
                    ndof-dimensional vector that represents the system's geometry. If provided (other than None),
                    it will be used to update the geometries of the system will  [ default: None ]
                * **pes_params["xmin"]** ( double ): minimal range of PES scan, matters only if coord_type = 0 [ default: -10.0 ]
                * **pes_params["xmax"]** ( double ): maximal range of PES scan, matters only if coord_type = 0 [ default:  10.0 ]
                * **pes_params["dx"]** ( double ): PES scan step size, matters only if coord_type = 0 [ default:  1.0 ]
                * **pes_params["tmin"]** ( int ): minimal range of t-PES scan, matters only if coord_type = 1 [ default: 0.0 ]
                * **pes_params["tmax"]** ( int ): maximal range of t-PES scan, matters only if coord_type = 1 [ default: 10.0 ]
                * **pes_params["dt"]** ( int ): t-PES scan step size, matters only if coord_type = 1 [ default:  1 ]
                * **pes_params["rep_tdse"]** ( 0 or 1 ): representation we are interested in
                    ( 0 - diabatic, 1 - adiabatic) [ required ]
                * **pes_params["rep_ham"]** ( 0 or 1 ): representation of the Hamiltonian returned 
                    by the `comp_model` function ( 0 - diabatic, 1 - adiabatic) [ required ]                    
                    
            plot_params_ ( dictionary ): determines what to print and how to do it, can contain the following keys:
            
                * **pes_params["which_ham_dia"]** ( list of 2-element lists ): which matrix elements of a 
                    diabatic Hamiltonian to plot [ default: empty ]
                    
                * **pes_params["which_ham_adi"]** ( list of 2-element lists ): which matrix elements of an 
                    adiabatic Hamiltonian to plot [ default: empty ]
                    
                * **pes_params["which_d1ham_dia"]** ( list of 3-element lists ): which matrix elements of a 
                    derivatives of a diabatic Hamiltonian w.r.t. which nuclear DOFs to plot, the format of each
                    entry is [idof, istate, jstate] [ default: empty ]

                * **pes_params["which_d1ham_adi"]** ( list of 3-element lists ): which matrix elements of a 
                    derivatives of an adiabatic Hamiltonian w.r.t. which nuclear DOFs to plot, the format of each
                    entry is [idof, istate, jstate] [ default: empty ]

                * **pes_params["which_dc1_dia"]** ( list of 3-element lists ): which matrix elements of a 
                    derivative couplings in a diabatic representation w.r.t. which nuclear DOFs to plot, the format
                    of each entry is [idof, istate, jstate] [ default: empty ]

                * **pes_params["which_dc1_adi"]** ( list of 3-element lists ): which matrix elements of a 
                    derivative couplings in an adiabatic representation w.r.t. which nuclear DOFs to plot, the format
                    of each entry is [idof, istate, jstate] [ default: empty ]

                * ** plot_params["colors"]** ( dictionary ) : defines the list of color definition (similar to the one found below)

                * ** plot_params["clrs_index"]** ( list of strings ): the mapping of a color definition to its "name" (similar to the one found below)

                    
    """
    
    pes_params = dict(pes_params_)
    pes_params_critical = [ "rep_tdse", "rep_ham" ] 
    pes_params_default = { "ndia":2, "nadi":2, "ndof":1,
                           "active_dof":0,
                           "coord_type":0,
                           "reference_coord":MATRIX(1, 1),
                           "coord_mapping":None,
                           "xmin":-10.0, "xmax":10.0, "dx":1.0,
                           "tmin":0, "tmax":10, "dt":1
                         }
    comn.check_input(pes_params, pes_params_default, pes_params_critical)

    active_dof = pes_params["active_dof"]
    coord_type = pes_params["coord_type"]
    reference_coord = pes_params["reference_coord"]
    coord_mapping = pes_params["coord_mapping"]
    xmax = pes_params["xmax"]
    xmin = pes_params["xmin"]
    dx = pes_params["dx"]
    tmax = pes_params["tmax"]
    tmin = pes_params["tmin"]
    dt = pes_params["dt"]
    ndia = pes_params["ndia"]
    nadi = pes_params["nadi"]
    ndof = pes_params["ndof"]
    rep_tdse = pes_params["rep_tdse"]
    rep_ham = pes_params["rep_ham"]
    



    colors_ = {}

    colors_.update({"11": "#8b1a0e"})  # red       
    colors_.update({"12": "#FF4500"})  # orangered 
    colors_.update({"13": "#B22222"})  # firebrick 
    colors_.update({"14": "#DC143C"})  # crimson   

    colors_.update({"21": "#5e9c36"})  # green
    colors_.update({"22": "#006400"})  # darkgreen  
    colors_.update({"23": "#228B22"})  # forestgreen
    colors_.update({"24": "#808000"})  # olive      

    colors_.update({"31": "#8A2BE2"})  # blueviolet
    colors_.update({"32": "#00008B"})  # darkblue  

    colors_.update({"41": "#2F4F4F"})  # darkslategray

    clrs_index_ = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]
    
    plot_params = dict(plot_params_)
    plot_params_critical = [  ] 
    plot_params_default = { "which_ham_dia":[],    "which_ham_adi":[],
                            "which_d1ham_dia":[],  "which_d1ham_adi":[],
                            "which_dc1_dia":[],    "which_dc1_adi":[],
                            "colors":colors_, "clrs_index":clrs_index_
                          }
    comn.check_input(plot_params, plot_params_default, plot_params_critical)
    which_ham_dia = plot_params["which_ham_dia"]
    which_ham_adi = plot_params["which_ham_adi"]
    which_d1ham_dia = plot_params["which_d1ham_dia"]
    which_d1ham_adi = plot_params["which_d1ham_adi"]
    which_dc1_dia = plot_params["which_dc1_dia"]
    which_dc1_adi = plot_params["which_dc1_adi"]
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    
    
        
    grid = []
    nsteps = 0
    if coord_type==0:
        nsteps = int((xmax - xmin) / dx) + 1
        for i in range(nsteps):
            grid.append(xmin + i * dx)
    elif coord_type==1:
        nsteps = int((tmax - tmin) / dt) + 1
        for i in range(nsteps):
            grid.append(tmin + i * dt)

        
            
    # ======= Hierarchy of Hamiltonians =======
    tol = 0.01
    ham = nHamiltonian(ndia, nadi, ndof)
    ham.init_all(2)
    ham.phase_corr_ovlp_tol = tol

    ham1 = [] 
    for tr in range(1):
        ham1.append( nHamiltonian(ndia, nadi, ndof) )        
        ham1[tr].init_all(2)
        ham1[tr].phase_corr_ovlp_tol = tol
        ham.add_child(ham1[tr])
        
    projectors = CMATRIXList()
    for tr in range(1):
        projectors.append(CMATRIX(nadi, nadi))
        projectors[tr].identity()   

        
    # Energies, forces, and couplings    
    ham_dia = []
    ham_adi = []
    d1ham_dia = []
    d1ham_adi = []
    dc1_dia = []
    dc1_adi = []
    
    for i in which_ham_dia:
        ham_dia.append([])
        
    for i in which_ham_adi:
        ham_adi.append([])
        
    for i in which_d1ham_dia:
        d1ham_dia.append([])
        
    for i in which_d1ham_adi:
        d1ham_adi.append([])
        
    for i in which_dc1_dia:
        dc1_dia.append([])
        
    for i in which_dc1_adi:
        dc1_adi.append([])
        
                    
    nsteps = len(grid)
    scan_path = []
    
    tid = Py2Cpp_int([0, 0])

    dyn_var = dyn_variables(ndia, nadi, ndof, 1)
    
    for step in range(nsteps):
        q = None
        if coord_mapping == None:
            q = MATRIX(reference_coord); 
            q.set(active_dof, 0, grid[step])
        else:
            q = coord_mapping(grid[step])
        scan_path.append( MATRIX(q) )

        model_params["timestep"] = step
            
        # Old way:
        # update_Hamiltonian_q( {"rep_tdse":rep_tdse, "rep_ham":rep_ham}, 
        #                      q, projectors, ham, comp_model, model_params)

        # New way
        dyn_var.set_q(q)
        update_Hamiltonian_variables( { "time_overlap_method":1, "ham_update_method":2, 
                                        "ham_transform_method":0, "nac_update_method":0,
                                        "hvib_update_method":0, "rep_tdse":rep_tdse, "rep_ham":rep_ham}, 
                                         dyn_var, ham, ham, comp_model, model_params, 0 ) 
        
        #if rep_tdse==0: 
            
        for it_indx, it in enumerate(which_ham_dia):                
            ham_dia[it_indx].append( ham.get_ham_dia(tid).get(it[0], it[1]).real )
                
        for it_indx, it in enumerate(which_d1ham_dia):                
            d1ham_dia[it_indx].append( ham.get_d1ham_dia(it[0],tid).get(it[1], it[2]).real )

        for it_indx, it in enumerate(which_dc1_dia):                
            dc1_dia[it_indx].append( ham.get_dc1_dia(it[0],tid).get(it[1], it[2]).real )
                
            
        #elif rep_tdse==1:   
            
        for it_indx, it in enumerate(which_ham_adi):                
            ham_adi[it_indx].append( ham.get_ham_adi(tid).get(it[0], it[1]).real )

        for it_indx, it in enumerate(which_d1ham_adi):                
            d1ham_adi[it_indx].append( ham.get_d1ham_adi(it[0],tid).get(it[1], it[2]).real )

        for it_indx, it in enumerate(which_dc1_adi):                
            dc1_adi[it_indx].append( ham.get_dc1_adi(it[0],tid).get(it[1], it[2]).real )

        

       
    plt.rc('axes', titlesize=38)      # fontsize of the axes title
    plt.rc('axes', labelsize=38)      # fontsize of the x and y labels
    plt.rc('legend', fontsize=36)     # legend fontsize
    plt.rc('xtick', labelsize=28)     # fontsize of the tick labels
    plt.rc('ytick', labelsize=28)     # fontsize of the tick labels

    plt.rc('figure.subplot', left=0.2)
    plt.rc('figure.subplot', right=0.95)
    plt.rc('figure.subplot', bottom=0.13)
    plt.rc('figure.subplot', top=0.88)


            
    #======== Now lets plot what we have computed ===========
    plt.figure(1, figsize=(36, 18)) # dpi=300, frameon=False)
    plt.subplot(1,2,1)
    plt.title('Diabatic energies')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Energy, a.u.')    
    for it_indx, it in enumerate(which_ham_dia):
        plt.plot(grid, ham_dia[it_indx], label='$H_{%i, %i}$' %(it[0], it[1]), linewidth=5, color = colors[ clrs_index[it_indx] ]) 
    plt.legend()
    
    plt.subplot(1,2,2)
    plt.title('Adiabatic energies')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Energy, a.u.')    
    for it_indx, it in enumerate(which_ham_adi):
        plt.plot(grid, ham_adi[it_indx], label='$H_{%i, %i}$' % (it[0], it[1]), linewidth=5, color = colors[ clrs_index[it_indx] ]) 
    plt.legend()
    plt.show()
    plt.close()

    plt.figure(2, figsize=(36, 18)) # dpi=300, frameon=False)
    plt.subplot(1,2,1)
    plt.title('Derivatives of diabatic energies')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Derivative energy, a.u.')    
    for it_indx, it in enumerate(which_d1ham_dia):
        plt.plot(grid, d1ham_dia[it_indx], label='$dH_{%i, %i} / dR_{%i}$' % (it[1],it[2], it[0]), linewidth=5, color = colors[ clrs_index[it_indx] ]) 
    plt.legend()
    
    plt.subplot(1,2,2)
    plt.title('Derivatives of adiabatic energies')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Derivative energy, a.u.')    
    for it_indx, it in enumerate(which_d1ham_adi):
        plt.plot(grid, d1ham_adi[it_indx], label='$dH_{%i, %i} / dR_{%i}$' % (it[1],it[2], it[0]), linewidth=5, color = colors[ clrs_index[it_indx] ]) 
    plt.legend()
    plt.show()
    plt.close()
    
    
    plt.figure(3, figsize=(36, 18)) # dpi=300, frameon=False)
    plt.subplot(1,2,1)
    plt.title('Derivatives couplings, diabatic rep.')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Derivative coupling, a.u.')    
    for it_indx, it in enumerate(which_dc1_dia):
        plt.plot(grid, dc1_dia[it_indx], label='$< \psi_{%i} | dR_{%i} | \psi_{%i} >$' % (it[1], it[0], it[2]), linewidth=5, color = colors[ clrs_index[it_indx] ]) 
    plt.legend()
    
    plt.subplot(1,2,2)
    plt.title('Derivatives couplings, adiabatic rep.')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Derivative energy, a.u.')    
    for it_indx, it in enumerate(which_dc1_adi):
        plt.plot(grid, dc1_adi[it_indx], label='$< \psi_{%i} | dR_{%i} | \psi_{%i} >$' % (it[1], it[0], it[2]), linewidth=5, color = colors[ clrs_index[it_indx] ]) 
    plt.legend()    
    

    plt.show()
    plt.close()

    return scan_path



def plot_surfaces(_compute_model, _param_sets, states_of_interest, xmin, xmax, dx, plot_params,\
                 _ndof=1, _active_dof=0, _all_coordinates=[0.0]):
    """
    Args:
        _compute_model ( PyObject ): the function that returns the class with Hamiltonian properties
        _param_sets ( list of lists of dictionaries ): parameters of the model Hamiltonian, many sets are possible (hense the list)

            For each set, the following keywords are required:
            * **nstates** ( int ): the dimensionality of the Hamiltonian

        states_of_interest ( list of ints ): indices of the states we want to plot 
        xmin ( double ): minimal value of the x axis used in the actual PES calculations  [a.u.]
        xmax ( double ): maximal value of the x axis used in the actual PES calculations  [a.u.]
        dx ( double ): step size of PES scan [a.u.]
        plot_params ( dictionary ): the parameters of plotting
            
            The dictionary can contain the following parameters:
            * **plot_params["colors"]** ( dictionary ) : defines the list of color definition (similar to the one found below)
            * **plot_params["clrs_index"]** ( list of strings ): the mapping of a color definition to its "name" (similar to the one found below)
            * **plot_params["xlim"]** ( list of 2 doubles ): the minimal and maximal values of the x axis in the plotted frame [ default: [-7.5, 15.0]]
            * **plot_params["ylim"]** ( list of 2 doubles ): the minimal and maximal values of the y axis in the plotted frame,
               for diabatic and adiabatic energies [default: [ -0.005, 0.025] ]
            * **plot_params["ylim2"]** ( list of 2 doubles ): the minimal and maximal values of the y axis in the plotted frame,
               for NACs [default: [-0.01, 0.01]]
            * **plot_params["do_show"]** (0 or 1) : 
               - 0: - don't show the plot in the display (e.g. when there is no tunneling or proper terminal)
               - 1: - do show it. Be careful - sometimes showing is not recommended because you'll get a bunch of windows which 
                 you'll need to close before the code would proceed to other parts. Turn it on if you use Jupyter notebooks. [ default: 1]
            * **plot_params["save_figures"]** (0 or 1): 0 - don't save the figures into separate files, 1 - do save them (e.g. when you can't or don't want
               to visualize them in terminal). The figures will be saved into a directory defined by the `prefix` variable. [default: 1]
            * **plot_params["prefix"]** (string): the name of the folder to which the figures will be saved (only if save_figures is set to 1). If the
               folder doesn't exist, it will be created [default: "out"]
            * **plot_params["dpi"]** (int): the quality of the figures to be saved [default: 300]
            * **plot_params["nac_idof"]** (int): when plotting derivative NACs, this variable defines the index of the nuclear DOF that defines such 
               a NAC, e.g. <\psi_i | d/dR_{nac_idof}| \psi_j > [default: 0]
            * **plot_params["plotting_option"]** (0 or 1): how to plot the results:
               - 0 : Plot diabatic and adiabatic surfaces separately, plot projections too [ default ]
               - 1 : Plot diabatic and adiabatic surfaces in one picture, using dashed lines for diabatic. Plot NACs in a separate panel
            * **plot_params["figsize"]** (list of ints) : define the size of the figure [ default: [36, 18]]
            * **plot_params["titlesize"]** (int): size of the title [default: 38]
            * **plot_params["labelsize"]** (int): size of the label [default: 38]
            * **plot_params["fontsize"]** (int): size of the fonts [default: 36]
            * **plot_params["xticksize"]** (int): size of the x-ticks [default: 28]
            * **plot_params["yticksize"]** (int): size of the y-ticks [default: 28]
            * **plot_params[""show_nac_abs]** (int): when plotting derivative NACs, this will also add plotting of the modulus of NACs:
               - 0 : do not show [default]
               - 1 : do show it
        
        _ndof ( int ): the dimensionality of the PES [ default: 1 ]
        _active_dof ( int ): the index of the DOF used to construct the PES [ default: 0 ]
        _all_coodinates ( list of doubles ): values of all coordinates, the one at the position of `_active_dof` will be disregarded, while all
            other will be fixed at the values provided
        
    """


    colors_ = {}

    colors_.update({"11": "#8b1a0e"})  # red       
    colors_.update({"12": "#FF4500"})  # orangered 
    colors_.update({"13": "#B22222"})  # firebrick 
    colors_.update({"14": "#DC143C"})  # crimson   

    colors_.update({"21": "#5e9c36"})  # green
    colors_.update({"22": "#006400"})  # darkgreen  
    colors_.update({"23": "#228B22"})  # forestgreen
    colors_.update({"24": "#808000"})  # olive      

    colors_.update({"31": "#8A2BE2"})  # blueviolet
    colors_.update({"32": "#00008B"})  # darkblue  

    colors_.update({"41": "#2F4F4F"})  # darkslategray

    clrs_index_ = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]



    # Parameters and dimensions
    critical_params = [  ] 
    default_params = {  "colors":colors_, "clrs_index":clrs_index_,
                        "xlim":[-7.5, 15], "ylim":[-0.005, 0.025], "ylim2":[-0.01, 0.01], "do_show":1,
                        "save_figures":1, "prefix":"out", "dpi":300, "nac_idof":0,
                        "plotting_option":0, "figsize":[36, 18], "titlesize":38, "labelsize":38,
                        "fontsize": 36, "xticksize":28, "yticksize":28, "show_nac_abs":0,
                     }
    comn.check_input(plot_params, default_params, critical_params)
        
    colors = plot_params["colors"]
    clrs_index = plot_params["clrs_index"]
    xlim = plot_params["xlim"]
    ylim = plot_params["ylim"] # range of diabatic and adiabatic energies
    ylim2 = plot_params["ylim2"] # range of NACs
    do_show = plot_params["do_show"]
    save_figures = plot_params["save_figures"]
    prefix = plot_params["prefix"]
    dpi_value = plot_params["dpi"]
    plotting_option = plot_params["plotting_option"]  
    nac_idof = plot_params["nac_idof"]  # index of the DOF for which we plot the derivative coupling
    _figsize = plot_params["figsize"]
    _titlesize = plot_params["titlesize"]
    _labelsize = plot_params["labelsize"]
    _fontsize = plot_params["fontsize"]
    _xticksize = plot_params["xticksize"]
    _yticksize = plot_params["yticksize"]
    show_nac_abs = plot_params["show_nac_abs"]

    fig_size = (_figsize[0], _figsize[1])

    
    X = []
    nsteps = int((xmax - xmin) / dx) + 1

    for i in range(nsteps):
        X.append(xmin + i * dx)
    
    
    plt.rc('axes', titlesize=_titlesize)      # fontsize of the axes title
    plt.rc('axes', labelsize=_labelsize)      # fontsize of the x and y labels
    plt.rc('legend', fontsize=_fontsize)      # legend fontsize
    plt.rc('xtick', labelsize=_xticksize)     # fontsize of the tick labels
    plt.rc('ytick', labelsize=_yticksize)     # fontsize of the tick labels

    plt.rc('figure.subplot', left=0.2)
    plt.rc('figure.subplot', right=0.95)
    plt.rc('figure.subplot', bottom=0.13)
    plt.rc('figure.subplot', top=0.88)
    

    sz = len(_param_sets)
    for iset in range(sz):
               
        # The "nstates" field must be specified
        comn.check_input(_param_sets[iset], {}, ["nstates"])
        n = _param_sets[iset]["nstates"]
        nstates = n
        
        ham = nHamiltonian(nstates, nstates, _ndof) # ndia, nadi, nnucl
        ham.init_all(2)

        
        hdia, hadi, nac, nac_abs  = [], [], [], []
        uij = []              # projecitions of the MOs onto elementary basis
    
        for k1 in range(nstates):
            hadi.append([])
            hdia.append([])
            nac_k1 = []
            nac_abs_k1 = []
            uij_k1 = []
            for k2 in range(nstates):
                uij_k1.append([])
                nac_k1.append([])
                nac_abs_k1.append([])
            uij.append(uij_k1)
            nac.append(nac_k1)
            nac_abs.append(nac_abs_k1)
                
    
        for i in range(nsteps):

            scan_coord = MATRIX(_ndof, 1);
            for j in range(_ndof):
                scan_coord.set(j, 0, _all_coordinates[j])                 
            scan_coord.set(_active_dof, 0, X[i])  
        
            # Diabatic properties
            ham.compute_diabatic(_compute_model, scan_coord, _param_sets[iset])
            
            # Adiabatic properties
            ham.compute_adiabatic(1);       
        
            U = ham.get_basis_transform()
            #P = U * U.H()  # population matrix
        
            for k1 in range(nstates):
                hadi[k1].append(ham.get_ham_adi().get(k1, k1).real)
                hdia[k1].append(ham.get_ham_dia().get(k1, k1).real)
            
                for k2 in range(nstates):
                    uij[k1][k2].append(U.get(k1,k2).real**2 + U.get(k1,k2).imag**2)
                    nac_k1_k2 = ham.get_dc1_adi(0).get(k1, k2).real
                    nac[k1][k2].append(nac_k1_k2)
                    nac_abs[k1][k2].append( abs(nac_k1_k2) )


        if plotting_option==0:  # Plot diabatic and adiabatic surfaces separately, plot projections too
                    
            plt.figure(2*iset, figsize=fig_size ) # dpi=300, frameon=False)
                
            plt.subplot(1, 2, 1)    
            plt.ylim(ylim[0], ylim[1])
            plt.xlim(xlim[0], xlim[1])
    
            #plt.title('Params set %i: Ham_dia' % (iset) )
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('Energy, a.u.')
            for k1 in states_of_interest:
                plt.plot(X, hdia[k1], label='$H_{%i%i}$' % (k1,k1), linewidth=7, color = colors[clrs_index[k1]])     
            plt.legend()    
    
            plt.subplot(1, 2, 2)
            plt.ylim(ylim[0], ylim[1])
            plt.xlim(xlim[0], xlim[1])
            #plt.title('Params set %i: Ham_adi' % (iset))
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('Energy, a.u.')
            for k1 in states_of_interest:
                plt.plot(X, hadi[k1], label='$E_{%i}$' % (k1), linewidth=7, color = colors[clrs_index[k1]])     
            plt.legend()    

            if save_figures:
                if not os.path.exists(F"{prefix}"):
                    os.system(F"mkdir {prefix}")
                plt.savefig(F"{prefix}/Ham_dia_E_adi_set_{iset}.png", dpi=dpi_value)            
            
                  
            plt.figure(2*iset+1, figsize=fig_size ) # dpi=300, frameon=False)            
            sz1 = len(states_of_interest)
    
            for k2 in states_of_interest:
                indx = states_of_interest.index(k2)

                plt.figure(2*iset+1+indx, figsize=(36, 18)) # dpi=300, frameon=False)            

                #plt.subplot(1, sz1, 1+indx)
                #plt.subplot(1, 1, 1+indx)
                #plt.title('Params set %i: Adi state %i' % (iset, k2) )
                plt.xlabel('Coordinate, a.u.')
                plt.ylabel('Projection')
    
                for k1 in range(nstates):
                    #plt.plot(X, uij[k1][k2], label='$dia_{%i} | adi_{%i}$' % (k1, k2), linewidth=7, color = colors[clrs_index[k1]])         
                    plt.plot(X, uij[k1][k2], label='$< \psi^{dia}_{%i} | \psi^{adi}_{%i} >$' % (k1, k2), linewidth=7, color = colors[clrs_index[k1]])         
                plt.legend()               

            if save_figures:
                if not os.path.exists(F"{prefix}"):
                    os.system(F"mkdir {prefix}")
                plt.savefig(F"{prefix}/projections_set_{iset}.png", dpi=dpi_value)            
    
            if do_show:
                plt.show()
        
            plt.close()          

        elif plotting_option==1:  # Plot diabatic and adiabatic surfaces in one picture, using dashed lines for diabatic
                                  # Plot NACs in a separate panel

            plt.figure(iset, figsize=fig_size) # dpi=300, frameon=False)

            plt.subplot(1, 2, 1)
            plt.ylim(ylim[0], ylim[1])
            plt.xlim(xlim[0], xlim[1])

            #plt.title('Params set %i: Ham_dia' % (iset) )
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('Energy, a.u.')

            for k1 in states_of_interest:
                plt.plot(X, hdia[k1], label='$H_{%i%i}$' % (k1,k1), linewidth=7, ls="--", color = colors[clrs_index[k1]])
                plt.plot(X, hadi[k1], label='$E_{%i}$' % (k1), linewidth=7, color = colors[clrs_index[k1]])
            plt.legend()


            plt.subplot(1, 2, 2)
            plt.ylim(ylim2[0], ylim2[1])
            plt.xlim(xlim[0], xlim[1])
            #plt.title('Params set %i: Ham_adi' % (iset))
            plt.xlabel('Coordinate, a.u.')
            plt.ylabel('NAC, a.u.')
            cnt = 0
            for k1 in states_of_interest:
                for k2 in states_of_interest:
                    if k2>k1:                        
                        plt.plot(X, nac[k1][k2], label='$NAC_{%i%i}$' % (k1,k2), linewidth=7, color = colors[clrs_index[cnt]])
                        if show_nac_abs:
                            plt.plot(X, nac_abs[k1][k2], label='$NAC_{%i%i}$' % (k1,k2), linewidth=7, ls="--", color = colors[clrs_index[cnt+1]])
                        cnt = cnt + 1
            plt.legend()


            if save_figures:
                if not os.path.exists(F"{prefix}"):
                    os.system(F"mkdir {prefix}")
                plt.savefig(F"{prefix}/Ham_dia_E_adi_NAC_set_{iset}.png", dpi=dpi_value)

            if do_show:
                plt.show()

            plt.close()
