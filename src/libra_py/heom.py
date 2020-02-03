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
.. module:: heom
   :platform: Unix, Windows
   :synopsis: This module implements a wrapper function for doing HEOM dynamics
       The code is a translation/refactoring of the Fortran code or Amber Jain & Joe Subotnik
       https://github.com/subotnikgroup/HEOM_Amber

       List of functions:
           * run_dynamics(dyn_params, Ham, rho_init)

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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


import util.libutil as comn
from . import units
from . import data_outs
from . import dynamics_io
from . import dynamics_hdf5


#===================== HEOM output ====================


def init_heom_data(saver, hdf5_output_level, _nsteps, _nquant):

    if hdf5_output_level>=1:
        # Time axis (integer steps)
        saver.add_dataset("timestep", (_nsteps,) , "I")  

        # Time axis
        saver.add_dataset("time", (_nsteps,) , "R")  
        

    if hdf5_output_level>=3:
        # System's density matrix
        saver.add_dataset("denmat", (_nsteps, _nquant, _nquant), "C") 



def save_heom_data(step, saver, params, denmat):
    
    dt = params["dt"]    
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level>=1:
        # Timestep 
        saver.save_scalar(step, "timestep", step) 

        # Actual time
        saver.save_scalar(step, "time", step * dt)        

    if hdf5_output_level>=3:
        # Average adiabatic density matrices
        saver.save_matrix(step, "denmat", denmat) 




def run_dynamics(dyn_params, Ham, rho_init):

    """
    
    Args: 
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        _p ( MATRIX(nnucl, ntraj) ): momenta of the "classical" particles [units: a.u. of momenta]
        _iM ( MATRIX(nnucl, 1) ): masses of classical particles [units: a.u.^-1]
        _Cdia ( CMATRIX(ndia, ntraj) ): amplitudes of the diabatic basis states
        _Cadi ( CMATRIX(nadi, ntraj) ): amplitudes of the adiabatic basis states
        _projectors ( list of CMATRIX(nadi, nadi) ): cumulative phase correction and state tracking matrices
        _states ( intList, or list of ntraj ints ): the quantum state of each trajectory

        dyn_params ( dictionary ): parameters controlling the execution of the dynamics
            Can contain:
      
            * **dyn_params["rep_tdse"]** ( int ): selects the representation in which 
                nuclear/electronic (Ehrenfest core) dynamics is executed. The representation 
                used to integrate TD-SE

                - 0: diabatic representation
                - 1: adiabatic representation [ default ]

    """


    params = dict(dyn_params)   

    # Parameters and dimensions
    critical_params = [  ] 
    default_params = { "KK":0, "LL":10,
                       "verbosity":-1,
                       "gamma": 1.0/(0.1 * units.ps2au),
                       "eta": 2.0 * 50.0 * units.inv_cm2Ha,
                       "temperature": 300.0,
                       "tolerance":1e-6,
                       "filter_after_steps":10,
                       "dt":0.1*units.fs2au, "nsteps":10, "progress_frequency":0.1,
                       "properties_to_save": [ "timestep", "time", "denmat"],
                       "use_compression":0, "compression_level":[0,0,0],
                       "hdf5_output_level":0, "prefix":"out", 
                       "txt_output_level":0, "mem_output_level":3,

                     }

    comn.check_input(params, default_params, critical_params)

    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"]*nsteps)    
    

    #============= System ======================
    params.update({"Ham" : Ham})
    nquant = Ham.num_of_cols

    #============== HEOM topology ==============


    KK = dyn_params["KK"]    
    LL = dyn_params["LL"] 


    nn_tot = compute_nn_tot(nquant, KK, LL)
    
    nn = allocate_3D(nquant+1, KK+1, nn_tot+1)
    map_nplus = allocate_3D(nquant+1, KK+1, nn_tot+1)
    map_nneg = allocate_3D(nquant+1, KK+1, nn_tot+1)
    zero = allocate_1D(nn_tot+1)
    map_sum = allocate_1D(LL+1)

    compute_nn(nquant, KK, LL, map_sum, nn);
    compute_map(nquant, KK, LL, nn, map_nplus, map_nneg);

    if params["verbosity"]>=0:
        print(F"nn_tot = {nn_tot}")


    params.update( { "nn":nn, "zero":zero, "map_nplus":map_nplus, "map_nneg":map_nneg } )


    #============ Bath update =====================
    gamma_matsubara = doubleList()
    c_matsubara = complexList()

    setup_bath(params, gamma_matsubara, c_matsubara)
    params.update({ "gamma_matsubara": gamma_matsubara, "c_matsubara":c_matsubara  } )

    if params["verbosity"]>=1:
        for k in range(KK+1):
            print(F" k = {k} gamma_matsubara[{k}] = {gamma_matsubara[k]}  c_matsubara[{k}] = {c_matsubara[k]}")

    #============= Initialization ============

    rho = CMATRIX((nn_tot+1)*nquant, nquant)  # all rho matrices stacked on top of each other
    rho_unpacked = CMATRIXList()

    for n in range(nn_tot+1):
        rho_unpacked.append( CMATRIX(nquant, nquant))


    # Initial conditions
    x_ = Py2Cpp_int(list(range(nquant, 2*nquant))) 
    y_ = Py2Cpp_int(list(range(nquant))) 
    push_submatrix(rho, rho_init, x_, y_)




    #================ Create savers ==================    
    prefix = params["prefix"]

    # Create an output directory, if not present    
    if not os.path.isdir(prefix):
        os.mkdir(prefix)

    properties_to_save = params["properties_to_save"]


    #====== HDF5 ========
    hdf5_saver = None
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level > 0:                
        hdf5_saver = dynamics_hdf5.hdf5_saver(F"{prefix}/data.hdf", properties_to_save) 
        hdf5_saver.set_compression_level(params["use_compression"], params["compression_level"])

        init_heom_data(hdf5_saver, hdf5_output_level, params["nsteps"], nquant)


    #====== TXT ========
    txt_saver = None
    if params["txt_output_level"] > 0:
        pass
    
    #====== MEM =========
    mem_saver = None
    mem_output_level = params["mem_output_level"]

    if mem_output_level > 0:
        mem_saver =  dynamics_hdf5.mem_saver(properties_to_save)
        init_heom_data(mem_saver, mem_output_level, params["nsteps"], nquant)

                         
    savers = {"hdf5_saver":hdf5_saver, "txt_saver":txt_saver, "mem_saver":mem_saver }
    



    unpack_rho(rho_unpacked, rho)


    #============== Propagation =============

    for step in range(params["nsteps"]):

        #================ Saving the data ==================
        if step%print_freq==0:
            print(F" step= {step}")
            
        # Save properties
        if savers["hdf5_saver"] != None:            
            save_heom_data(step, savers["hdf5_saver"], params, rho_unpacked[1])
            
        if savers["txt_saver"] != None:            
            pass
            
        if savers["mem_saver"] != None:            
            prms = dict(params)
            prms["hdf5_output_level"] = prms["mem_output_level"]
            save_heom_data(step, savers["mem_saver"], prms, rho_unpacked[1])



        if step % params["filter_after_steps"] == 1:
            unpack_rho(rho_unpacked, rho)
            params["zero"] = filter(rho_unpacked, params["tolerance"]);
            pack_rho(rho_unpacked, rho)

        rho = RK4(rho, params["dt"], compute_heom_derivatives, params)


    
    if savers["mem_saver"] != None:        
        savers["mem_saver"].save_data( F"{prefix}/mem_data.hdf", properties_to_save, "w")
        return mem_saver




