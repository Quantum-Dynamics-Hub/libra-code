#*********************************************************************************                     
#* Copyright (C) 2020 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: compute
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing exact on-the-grid dynamics
       List of functions:
           * init_wfc(params, _potential, model_params )
           * run_dynamics(wfc, params, model_params, savers)
           * run_relaxation(_params, _potential, model_params)
           * plot_mem(res, _params, model_params, plot_params)
           * plot_hdf5(plot_params)

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
import time

import h5py
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
import libra_py.units as units
import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers

from . import save



def init_wfc(params, _potential, model_params ):
    """
    This is a generic and a bit excessive initialization procedure - some of the parameters
    may not be needed for some of the methods, but it allows us to have a simple structure of the program
    
    istate = [rep, index_of_the_state]
    
    """
            
    critical_params = []
    default_params = { "nsteps":200, "dt":10.0, "progress_frequency":0.1,
                       "rmin":[-15.0], "rmax":[15.0], "dx":[0.1], "nstates":2,
                       "x0":[0.0], "p0":[0.0], "istate":[1,0], "masses":[2000.0], "k":[0.001]
                      }
    comn.check_input(params, default_params, critical_params)
    
    # Grid properties 
    dx = Py2Cpp_double(params["dx"])
    rmin = Py2Cpp_double(params["rmin"])
    rmax = Py2Cpp_double(params["rmax"])
    nstates = params["nstates"]  

    # Dynamical properties 
    nsteps = params["nsteps"]
    dt = params["dt"]    
    
    # Properties of the initial wavefunction
    istate = params["istate"]          
    k = params["k"]
    masses = Py2Cpp_double(params["masses"])
    

    x0 = Py2Cpp_double(params["x0"])
    p0 = Py2Cpp_double(params["p0"])
    ndof = len(x0)
    
    dx0_tmp = []
    for idof in range(ndof):
        sigmax2 = 0.5*math.sqrt(1.0/(k[idof]*masses[idof]))        
        dx0_tmp.append(math.sqrt(sigmax2))
    dx0 = Py2Cpp_double(dx0_tmp)    
        
        
    
    # Here we initialize the grid and wavefunction
    wfc = Wfcgrid2( rmin, rmax,  dx, nstates)
    
    wfc.direct_allocate_tmp_vars(0)    # last time-step wavefunction, diabatic rep
    wfc.direct_allocate_tmp_vars(1)    # last time-step wavefunction, adiabatic rep

    
    wfc.update_Hamiltonian(_potential, model_params, 0)  # update Hamiltonian: diabatic -
                                                        # need to compute diabatic-to-adiabatic transform
    wfc.update_Hamiltonian(_potential, model_params, 1)  # update Hamiltonian: adiabatic, NACs
        
    wfc.update_propagator_H(0.5*dt)                     # copute the dia-to-adi transform + exp(-i* V *dt/2)
    
    wfc.add_wfc_Gau(x0, p0, dx0, istate[1], 1.0+0.0j, istate[0])   # Add to the diabatic or adiabatic state
        
    if istate[0]==0:
        # If we added a diabatic state:
        wfc.update_adiabatic()        # update adiabatic wavefunction    
            
    elif istate[0]==1:
        # If we added an adiabatic state:
        wfc.update_diabatic()        # update diabatic wavefunction            
    
    
    wfc.update_propagator_K(dt, masses) # update reci space propagator in diabatic rep, exp(-iTdt) 
    wfc.update_reciprocal(1)     # update reci of adiabatic function        
    wfc.update_reciprocal(0)     # update reci of diabatic function
        
        
    print( "Norm (dia) = ", wfc.norm(0) )
    print( "Norm (adi) = ", wfc.norm(1) )
    print( "Ekin (dia) = ", wfc.e_kin(masses, 0) )
    print( "Ekin (adi) = ", wfc.e_kin(masses, 1) )
    print( "Epot (dia) = ", wfc.e_pot(0) )
    print( "Epot (adi) = ", wfc.e_pot(1) )

    
    return wfc





def run_dynamics(wfc, _params, model_params, savers):
    """
    Args:
        wfc ( Wfcgrid2 object ): object that stores the wavefunction and all related properties
        params ( dict ): parameters controlling the execution of the dynamics
            Can contain:

            * **params["nsteps"]** ( int ): the number of integration steps to compute [ default: 200 ]

            * **params["dt"]** ( float ): integration time steps [ a.u. of time; default: 10.0 ]
 
            * **params["progress_frequency"]** ( float ): the fraction of the total trajectory length 
               at which the progress status is printed out - just to monitor how the calculation is going [ default: 0.1 ]

            * **params["masses"]** (list of floats): the length of this list is equal to the number of nuclear DOFs 
               the list contains the masses of all DOFs [ units: atomic units; default: [2000.0] ]

            * **params["integrator"]** (string): the name of the integrator to use. Options:

                - "SOFT": Split-Operator Fourier Transform method of Kosloff and Kosloff [ default ]
                - "direct_dia": numeric, finite difference in diabatic representation 
                - "direct_adi": numeric, finite difference in adiabatic representation
                - "Colbert_Miller_dia": Colbert-Miller representation of kinetic energy, finite difference in diabatic representation
                - "Colbert_Miller_adi": Colbert-Miller representation of kinetic energy, finite difference in adiabatic representation
                - "Colbert_Miller_SOFT": Colbert-Miller representation of kinetic energy, using SOFT approach although the kinetic energy is 
                  computed in real space

            * **params["wfc_prefix"]** (string): the name of the folder to which we want save the wavefunction snapshots [ default: "wfc"]
 
            * **params["snap_freq"]** (int): the number of steps between each printout of the wfc snapshot [ default: 1 ]

            * **params["wfcr_params"]** ( list of 3 ints): each int is 0 or 1; 1 - meaning to print the corresponding property and 0 - not
                the flags are as follows: first - real part of the wfc, second - imaginary part of wfc, third - density of wfc. Indexing 
                starts at zero. This is all for the real-space wavefunction. [ default: [0] ]

            * **params["wfcr_rep"]** (int): representation for which the real-space wfc is printed out, 0 - diabatic, 1 - adiabatic [ default: 0 ]

            * **params["wfcr_states"]** (list of ints): indices of the states for which to print the corresponding real-space wavefunctions

            * **params["wfck_params"]** ( list of 3 ints): each int is 0 or 1; 1 - meaning to print the corresponding property and 0 - not
                the flags are as follows: first - real part of the wfc, second - imaginary part of wfc, third - density of wfc. Indexing
                starts at zero. This is all for the reciprocal-space wavefunction. [ default: [0] ]

            * **params["wfck_rep"]** (int): representation for which the reciprocal-space wfc is printed out, 0 - diabatic, 1 - adiabatic [ default: 0 ]

            * **params["wfck_states"]** (list of ints): indices of the states for which to print the corresponding reciprocal-space wavefunctions

    """

    params = dict(_params)
    
    integrators_map = {"SOFT": 0,
                       "direct_dia": 1,
                       "direct_adi": 2,
                       "Colbert_Miller_dia":3,
                       "Colbert_Miller_adi":4,
                       "Colbert_Miller_SOFT":5
                      }


    critical_params = []
    default_params = { "nsteps":200, "dt":10.0, "progress_frequency":0.1, "masses":[2000.0], "integrator":"SOFT", 
                       "wfc_prefix":"wfc", "snap_freq":1, "wfcr_params":[1, 1, 1], "wfcr_rep":0, "wfcr_states":[0],
                       "wfck_params":[1, 1, 1], "wfck_rep":0, "wfck_states":[0]
                      }
    comn.check_input(params, default_params, critical_params)



    integrator_id = integrators_map[ params["integrator"] ]
        
    nsteps = params["nsteps"]      
    print_freq = int(params["progress_frequency"]*nsteps)    
    dt = params["dt"]    
    masses = Py2Cpp_double(params["masses"])    
    ndof = len(masses)
    wfc_prefix = params["wfc_prefix"]    
    snap_freq = params["snap_freq"]

    wfcr_params = params["wfcr_params"]
    wfcr_rep = params["wfcr_rep"]            
    wfcr_states = Py2Cpp_int(params["wfcr_states"])

    wfck_params = params["wfck_params"]
    wfck_rep = params["wfck_rep"]
    wfck_states = Py2Cpp_int(params["wfck_states"])


    #================ Special setups ===================
    expT, expV = None, None
    if integrator_id == 5:
       
        T = wfc.operator_T(Py2Cpp_int([1]), masses, 1.0+0.0j) 
    
        expT = CMATRIX(wfc.Npts, wfc.Npts)
        exp_matrix(expT, T, -dt*1.0j, 0)
    
        expV = CMATRIXList()
        for ipt in range(wfc.Npts):
            expv = CMATRIX(wfc.nstates, wfc.nstates)
            exp_matrix(expv, wfc.Hdia[ipt], -dt*0.5j, 0);        
            expV.append(CMATRIX(expv))
        
    if not os.path.exists(wfc_prefix):
        os.system(F"mkdir {wfc_prefix}")
        
    #=============== Propagation ==========================
        
    for step in range(nsteps):
        
        #================ Saving the data ==================
        if step%print_freq==0:
            print(F" step= {step}")
            
        # Save properties
        if savers["hdf5_saver"] != None:            
            save.save_data_hdf5(step, wfc, savers["hdf5_saver"], params)
            
        if savers["txt_saver"] != None:            
            pass
            #save_data_txt(step, wfc, savers["txt_saver"], params)
            
        if savers["mem_saver"] != None:            
            prms = dict(params)
            prms["hdf5_output_level"] = prms["mem_output_level"]
            save.save_data_hdf5(step, wfc, savers["mem_saver"], prms)
    
        #================ Integration ==================
    
        if integrator_id==0:  # SOFT            
            wfc.SOFT_propagate()      # evolve the diabatic wfc 
            
        elif integrator_id == 1: # direct_dia
            if step==0:            
                wfc.direct_propagate_dia1(dt, masses)                           
            else:
                wfc.direct_propagate_dia2(dt, masses)                          
                
        elif integrator_id == 2: # direct_adi
            if step==0:            
                wfc.direct_propagate_adi1(dt, masses)                           
            else:
                wfc.direct_propagate_adi2(dt, masses)                                          
                
        elif integrator_id == 3: # Colbert-Miller_dia        
            if step==0:            
                wfc.Colbert_Miller_propagate_dia1(dt, masses)                           
            else:
                wfc.Colbert_Miller_propagate_dia2(dt, masses)                          

        elif integrator_id == 4: # Colbert-Miller_adi        
            if step==0:            
                wfc.Colbert_Miller_propagate_adi1(dt, masses)                           
            else:
                wfc.Colbert_Miller_propagate_adi2(dt, masses)                                          

        elif integrator_id == 5: # Colbert-Miller_SOFT        
            wfc.Colbert_Miller_SOFT(expT, expV, 0)    
            
            
        #============= Update other variables ================
        if integrator_id in [0, 1, 3, 5]:
            wfc.update_adiabatic()  
            
        if integrator_id in [2, 4]:
            wfc.update_diabatic()              
        
        wfc.update_reciprocal(0)  # update reci of diabatic function
        wfc.update_reciprocal(1)  # update reci of adiabatic function
        


        #============ Printing wfc into files ==================
        # Print out the wavefunctions and probability densities
        if step%snap_freq==0:
            snap = int(step/snap_freq)
            if ndof==1:
                a, b, c = wfcr_params[0], wfcr_params[1], wfcr_params[2]
                wfc.print_wfc_1D(F"{wfc_prefix}/wfcr_snap_{snap}", wfcr_rep, wfcr_states, a,b,c)        

                a, b, c = wfck_params[0], wfck_params[1], wfck_params[2]
                wfc.print_reci_wfc_1D(F"{wfc_prefix}/wfck_snap_{snap}", wfck_rep, wfck_states, a,b,c)   

            elif ndof==2:
                for st in wfcr_states:
                    a, b, c, st = wfcr_params[0], wfcr_params[1], wfcr_params[2], wfcr_states[0]
                    wfc.print_wfc_2D(F"{wfc_prefix}/wfcr_snap_{snap}_state_{st}", wfcr_rep, st, a,b,c)

                for st in wfck_states:
                    a, b, c, st = wfck_params[0], wfck_params[1], wfck_params[2], wfck_states[0]
                    wfc.print_reci_wfc_2D(F"{wfc_prefix}/wfck_snap_{snap}_state_{st}", wfck_rep, st, a,b,c)           


def run_relaxation(_params, _potential, model_params):
    """

    # For diabatic populations in boxes
    _params["customs_pops"] = [ [0, [-100], [-10] ], 
                                [0, [-10], [10] ],
                                [0, [10], [100] ]
                              ]

    # For adiabatic populations in boxes
    _params["customs_pops"] = [ [1, [-100], [-10] ], 
                                [1, [-10], [10] ],
                                [1, [10], [100] ]
                              ]

    """
    
        
    params = dict(_params)                
    nstates = len(model_params["E_n"])
    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"]*nsteps)    
    params.update({"nstates":nstates})


    ncustom_pops = 0
    if "custom_pops" not in params.keys():
        ncustom_pops = 0
    else:
        if params["custom_pops"]!=None:
            ncustom_pops = len(params["custom_pops"])
            print(F"ncustom_pops = {ncustom_pops}")
            
    
    print("Run calculations with the dynamical parameters ", params)
    print("Run calculations with the model parameters ", model_params)    
    
        
    wfc = init_wfc(params, _potential, model_params )
    ndof = wfc.ndof
    ngrid = wfc.Npts
    
    
    #================ Create savers ==================    
    # Create an output directory, if not present  
    prefix = params["prefix"]
    
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
            
    # Simulation parameters                    
    f = open(F"{prefix}/_dyn_params.txt","w")
    f.write( str(params) );  f.close()
    
    f = open(F"{prefix}/_model_params.txt","w")
    f.write( str(model_params) );  f.close()    

    
    
    properties_to_save = params["properties_to_save"]
    
    #====== HDF5 ========
    hdf5_saver = None
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level > 0:                
        hdf5_saver = data_savers.hdf5_saver(F"{prefix}/data.hdf", properties_to_save) 
        hdf5_saver.set_compression_level(params["use_compression"], params["compression_level"])
        save.exact_init_hdf5(hdf5_saver, hdf5_output_level, nsteps, ndof, nstates, ngrid)
        if ncustom_pops > 0:
            save.exact_init_custom_hdf5(hdf5_saver, nsteps, ncustom_pops, nstates)  # boxed populations on adiabatic/diabatic states

    #====== TXT ========
    txt_saver = None
    if params["txt_output_level"] > 0:
        pass
    
    #====== MEM =========
    mem_saver = None
    mem_output_level = params["mem_output_level"]

    if mem_output_level > 0:
        mem_saver =  data_savers.mem_saver(properties_to_save)
        save.exact_init_hdf5(mem_saver, mem_output_level, nsteps, ndof, nstates, ngrid)
        if ncustom_pops > 0:
            save.exact_init_custom_hdf5(mem_saver, nsteps, ncustom_pops, nstates)  # boxed populations on adiabatic/diabatic states

                         
    savers = {"hdf5_saver":hdf5_saver, "txt_saver":txt_saver, "mem_saver":mem_saver }
    

    #==================== Dynamics ======================    
    start = time.time()                               
    
    res = run_dynamics(wfc, params, model_params, savers)                     
    
    end = time.time()    
    print(F"Calculation time = {end - start} seconds")
    
    
    if mem_saver != None:        
        mem_saver.save_data( F"{prefix}/mem_data.hdf", properties_to_save, "w")
        return mem_saver

        

