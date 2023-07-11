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
.. module:: savers
   :platform: Unix, Windows
   :synopsis: This module implements functions for initializing and storing data
       computed/produced during TSH/Ehrenfest/Adiabatic MD calculations. If you need to print out more 
       data, you shall add the corresponding info in these modules

       List of functions:

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
#import libra_py.units as units
#import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers


#===================== TSH calculations output ====================

def init_tsh_data(saver, output_level, _nsteps, _ntraj, _ndof, _nadi, _ndia):
    """
    saver - can be either hdf5_saver or mem_saver

    """
    #print(F"In init_tsh_data with:nsteps = {_nsteps}, ntraj = {_ntraj}, ndof = {_ndof}, nadi = {_nadi}, ndia = {_ndia} ")
    #print(F"output_level = {output_level}")

    if output_level>=1:

        # Time axis (integer steps)   
        if "timestep" in saver.keywords:  
            saver.add_dataset("timestep", (_nsteps,) , "I")  

        # Time axis
        if "time" in saver.keywords:
            saver.add_dataset("time", (_nsteps,) , "R")  
        
        # Average kinetic energy
        if "Ekin_ave" in saver.keywords:
            saver.add_dataset("Ekin_ave", (_nsteps,) , "R")  
        
        # Average potential energy
        if "Epot_ave" in saver.keywords: 
            saver.add_dataset("Epot_ave", (_nsteps,) , "R")  
        
        # Average total energy
        if "Etot_ave" in saver.keywords:
            saver.add_dataset("Etot_ave", (_nsteps,) , "R")  
        
        # Fluctuation of average kinetic energy
        if "dEkin_ave" in saver.keywords:
            saver.add_dataset("dEkin_ave", (_nsteps,) , "R")  
        
        # Fluctuation of average potential energy
        if "dEpot_ave" in saver.keywords:
            saver.add_dataset("dEpot_ave", (_nsteps,) , "R")  
        
        # Fluctuation of average total energy
        if "dEtot_ave" in saver.keywords:
            saver.add_dataset("dEtot_ave", (_nsteps,) , "R")  

        # Thermostat energy 
        if "Etherm" in saver.keywords:
            saver.add_dataset("Etherm", (_nsteps,) , "R")  

        # System + thermostat energy
        if "E_NHC" in saver.keywords:
            saver.add_dataset("E_NHC", (_nsteps,) , "R")  


    if output_level>=2:

        # Trajectory-resolved instantaneous adiabatic states
        if "states" in saver.keywords: # and "states" in saver.np_data.keys():
            saver.add_dataset("states", (_nsteps, _ntraj), "I") 

        # Average adiabatic SE populations 
        if "se_pop_adi" in saver.keywords: # and "states" in saver.np_data.keys():
            saver.add_dataset("se_pop_adi", (_nsteps, _nadi), "R") 

        # Average diabatic SE populations 
        if "se_pop_dia" in saver.keywords: # and "states" in saver.np_data.keys():
            saver.add_dataset("se_pop_dia", (_nsteps, _ndia), "R") 

        # Average adiabatic SH populations 
        if "sh_pop_adi" in saver.keywords: # and "states" in saver.np_data.keys():
            saver.add_dataset("sh_pop_adi", (_nsteps, _nadi), "R") 

        # Average diabatic SH populations 
        if "sh_pop_dia" in saver.keywords: # and "states" in saver.np_data.keys():
            saver.add_dataset("sh_pop_dia", (_nsteps, _ndia), "R") 



    if output_level>=3:

        # Average adiabatic SH populations (dynamically-consistent)
        if "SH_pop" in saver.keywords: # and "SH_pop" in saver.np_data.keys():
            saver.add_dataset("SH_pop", (_nsteps, _nadi, 1), "R") 

        # Average adiabatic SH populations (raw)
        if "SH_pop_raw" in saver.keywords: # and "SH_pop_raw" in saver.np_data.keys():
            saver.add_dataset("SH_pop_raw", (_nsteps, _nadi, 1), "R") 


        # Average adiabatic density matrices (dynamically-consistent)
        if "D_adi" in saver.keywords: # and "D_adi" in saver.np_data.keys():
            saver.add_dataset("D_adi", (_nsteps, _nadi, _nadi), "C") 

        # Average adiabatic density matrices (raw)
        if "D_adi_raw" in saver.keywords: # and "D_adi_raw" in saver.np_data.keys():
            saver.add_dataset("D_adi_raw", (_nsteps, _nadi, _nadi), "C") 


        # Average diabatic density matrices (dynamically-consistent)
        if "D_dia" in saver.keywords: # and "D_dia" in saver.np_data.keys():
            saver.add_dataset("D_dia", (_nsteps, _ndia, _ndia), "C") 

        # Average diabatic density matrices (raw)
        if "D_dia_raw" in saver.keywords: # and "D_dia_raw" in saver.np_data.keys():
            saver.add_dataset("D_dia_raw", (_nsteps, _ndia, _ndia), "C") 



        # Trajectory-resolved coordinates
        if "q" in saver.keywords: # and "q" in saver.np_data.keys():
            saver.add_dataset("q", (_nsteps, _ntraj, _ndof), "R") 

        # Trajectory-resolved momenta
        if "p" in saver.keywords: # and "p" in saver.np_data.keys():
            saver.add_dataset("p", (_nsteps, _ntraj, _ndof), "R") 

        # Trajectory-resolved active forces
        if "f" in saver.keywords: # and "f" in saver.np_data.keys():
            saver.add_dataset("f", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved adiabatic TD-SE amplitudes
        if "Cadi" in saver.keywords: # and "Cadi" in saver.np_data.keys():
            saver.add_dataset("Cadi", (_nsteps, _ntraj, _nadi), "C") 

        # Trajectory-resolved diabatic TD-SE amplitudes
        if "Cdia" in saver.keywords: # and "Cdia" in saver.np_data.keys():
            saver.add_dataset("Cdia", (_nsteps, _ntraj, _ndia), "C") 


    if output_level>=4:

        # Trajectory-resolved vibronic Hamiltoninans in the adiabatic representation
        if "hvib_adi" in saver.keywords: # and "hvib_adi" in saver.np_data.keys():
            saver.add_dataset("hvib_adi", (_nsteps, _ntraj, _nadi, _nadi), "C") 

        # Trajectory-resolved vibronic Hamiltoninans in the diabatic representation
        if "hvib_dia" in saver.keywords: # and "hvib_dia" in saver.np_data.keys():
            saver.add_dataset("hvib_dia", (_nsteps, _ntraj, _ndia, _ndia), "C") 

        # Trajectory-resolved time-overlaps of the adiabatic states
        if "St" in saver.keywords: # and "St" in saver.np_data.keys():
            saver.add_dataset("St", (_nsteps, _ntraj, _nadi, _nadi), "C") 

        # Trajectory-resolved diabatic-to-adiabatic transformation matrices 
        if "basis_transform" in saver.keywords: # and "basis_transform" in saver.np_data.keys():
            saver.add_dataset("basis_transform", (_nsteps, _ntraj, _ndia, _nadi), "C") 

        # Trajectory-resolved projector matrices (from the raw adiabatic to consistent adiabatic)
        if "projector" in saver.keywords: # and "projector" in saver.np_data.keys():
            saver.add_dataset("projector", (_nsteps, _ntraj, _nadi, _nadi), "C") 





def init_tsh_savers(params, model_params, nsteps, ntraj, nnucl, nadi, ndia):

    #================ Create savers ==================    
    if params["txt_output_level"] > 0 or params["hdf5_output_level"] > 0 or params["mem_output_level"] > 0:

        prefix = params["prefix"]

        # Create an output directory, if not present    
        if not os.path.isdir(prefix):
            os.mkdir(prefix)

        # Simulation parameters                    
        f = open(F"{prefix}/_dyn_params.txt","w")
        f.write( str(params) );  f.close()
    
        f = open(F"{prefix}/_model_params.txt","w")
        f.write( str(model_params) );  f.close()    


    if params["txt2_output_level"] > 0:

        prefix2 = params["prefix2"]

        # Create an output directory, if not present    
        if not os.path.isdir(prefix2):
            os.mkdir(prefix2)

        # Simulation parameters                    
        f = open(F"{prefix2}/_dyn_params.txt","w")
        f.write( str(params) );  f.close()
    
        f = open(F"{prefix2}/_model_params.txt","w")
        f.write( str(model_params) );  f.close()    


    properties_to_save = params["properties_to_save"]


    _savers = {"hdf5_saver":None, "txt_saver":None, "mem_saver":None , "txt2_saver":None }

    #====== HDF5 ========
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level > 0:                
        _savers["hdf5_saver"] = data_savers.hdf5_saver(F"{prefix}/data.hdf", properties_to_save) 
        _savers["hdf5_saver"].set_compression_level(params["use_compression"], params["compression_level"])
        init_tsh_data(_savers["hdf5_saver"], hdf5_output_level, nsteps, ntraj, nnucl, nadi, ndia)


    #====== TXT ========
    txt_output_level = params["txt_output_level"]

    if params["txt_output_level"] > 0:
        _savers["txt_saver"] = data_savers.mem_saver(properties_to_save) 
        #_savers["txt_saver"].set_compression_level(params["use_compression"], params["compression_level"])
        init_tsh_data(_savers["txt_saver"], txt_output_level, nsteps, ntraj, nnucl, nadi, ndia)


    #====== TXT2: No intermediate memory allocation ========
    txt2_output_level = params["txt2_output_level"]

    if params["txt2_output_level"] > 0:
        _savers["txt2_saver"] = data_savers.mem_saver(properties_to_save) 

        # Here, nsteps is set to 1 since in this type of saver, we only care about the current values,
        # not all the timesteps - that would be too consuming
        init_tsh_data(_savers["txt2_saver"], txt2_output_level, 1, ntraj, nnucl, nadi, ndia)   


    
    #====== MEM =========
    mem_output_level = params["mem_output_level"]

    if mem_output_level > 0:
        _savers["mem_saver"] =  data_savers.mem_saver(properties_to_save)
        init_tsh_data(_savers["mem_saver"], mem_output_level, nsteps, ntraj, nnucl, nadi, ndia)


    return _savers                         
    




def save_hdf5_1D(saver, i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type==0:
        t = i

    # Timestep 
    saver.save_scalar(t, "timestep", i) 

    # Actual time
    saver.save_scalar(t, "time", dt*i)  

    # Average kinetic energy
    saver.save_scalar(t, "Ekin_ave", Ekin)  

    # Average potential energy
    saver.save_scalar(t, "Epot_ave", Epot)  

    # Average total energy
    saver.save_scalar(t, "Etot_ave", Etot)  

    # Fluctuation of average kinetic energy
    saver.save_scalar(t, "dEkin_ave", dEkin)  

    # Fluctuation of average potential energy
    saver.save_scalar(t, "dEpot_ave", dEpot)  

    # Fluctuation average total energy
    saver.save_scalar(t, "dEtot_ave", dEtot)  

    # Thermostat energy 
    saver.save_scalar(t, "Etherm", Etherm)  

    # System + thermostat energy
    saver.save_scalar(t, "E_NHC", E_NHC)  



def save_hdf5_1D_new(saver, i, params, dyn_var, ham, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    dt = params["dt"]

    Ekin, Epot, Etot = 0.0, 0.0, 0.0
    dEkin, dEpot, dEtot, Etherm = 0.0, 0.0, 0.0, 0.0
    E_NHC = 0.0

    t = 0
    if txt_type==0:
        t = i

    # Timestep 
    saver.save_scalar(t, "timestep", i) 

    # Actual time
    saver.save_scalar(t, "time", dt*i)  

    # Average kinetic energy
    Ekin = dyn_var.compute_average_kinetic_energy()
    saver.save_scalar(t, "Ekin_ave", Ekin)  

    # Average potential energy
    Epot = average_potential_energy(params, dyn_var, ham)
    saver.save_scalar(t, "Epot_ave", Epot)  

    # Average total energy
    Etot = Ekin + Epot
    saver.save_scalar(t, "Etot_ave", Etot)  

    # Fluctuation of average kinetic energy
    saver.save_scalar(t, "dEkin_ave", dEkin)  

    # Fluctuation of average potential energy
    saver.save_scalar(t, "dEpot_ave", dEpot)  

    # Fluctuation average total energy
    saver.save_scalar(t, "dEtot_ave", dEtot)  

    # Thermostat energy 
    saver.save_scalar(t, "Etherm", Etherm)  

    # System + thermostat energy
    saver.save_scalar(t, "E_NHC", E_NHC)  




def save_hdf5_2D(saver, i, states, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type==0:
        t = i

    # Trajectory-resolved instantaneous adiabatic states
    # Format: saver.add_dataset("states", (_nsteps, _ntraj), "I")        
    ntraj = len(states)
    for itraj in range(ntraj):
        saver.save_multi_scalar(t, itraj, "states", states[itraj])




def save_hdf5_2D_new(saver, i, dyn_var, ham, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type==0:
        t = i


    if "states" in saver.keywords and "states" in saver.np_data.keys():
        # Trajectory-resolved instantaneous adiabatic states
        # Format: saver.add_dataset("states", (_nsteps, _ntraj), "I")        
        ntraj = dyn_var.ntraj
        for itraj in range(ntraj):
            saver.save_multi_scalar(t, itraj, "states", dyn_var.act_states[itraj])

    if "se_pop_dia" in saver.keywords and "se_pop_dia" in saver.np_data.keys():
        # Average diabatic SE populations 
        # Format: saver.add_dataset("se_pop_dia", (_nsteps, _ndia), "R") 
        pops_se0 = dyn_var.compute_average_se_pop(0)
        ndia = dyn_var.ndia
        for ist in range(ndia):
            saver.save_multi_scalar(t, ist, "se_pop_dia", pops_se0[ist])

    if "se_pop_adi" in saver.keywords and "se_pop_adi" in saver.np_data.keys():
        # Average adiabatic SE populations 
        # Format: saver.add_dataset("se_pop_adi", (_nsteps, _nadi), "R") 
        pops_se1 = dyn_var.compute_average_se_pop(1)
        nadi = dyn_var.nadi
        for ist in range(nadi):
            saver.save_multi_scalar(t, ist, "se_pop_adi", pops_se1[ist])
    
    if "sh_pop_adi" in saver.keywords and "sh_pop_adi" in saver.np_data.keys():    
        # Average adiabatic SH populations 
        # Format: saver.add_dataset("sh_pop_adi", (_nsteps, _nadi), "R") 
        pops_sh = dyn_var.compute_average_sh_pop()
        nadi = dyn_var.nadi
        for ist in range(nadi):
            saver.save_multi_scalar(t, ist, "sh_pop_adi", pops_sh[ist])
        



def save_hdf5_3D(saver, i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type==0:
        t = i

    # Average adiabatic SH populations (dynamically-consistent)
    # Format: saver.add_dataset("SH_pop", (_nsteps, _nadi, 1), "R") 
    if "SH_pop" in saver.keywords and "SH_pop" in saver.np_data.keys():
        saver.save_matrix(t, "SH_pop", pops) 

    # Average adiabatic SH populations (raw)
    # Format: saver.add_dataset("SH_pop_raw", (_nsteps, _nadi, 1), "R") 
    if "SH_pop_raw" in saver.keywords and "SH_pop_raw" in saver.np_data.keys():
        saver.save_matrix(t, "SH_pop_raw", pops_raw) 



    # Average adiabatic density matrices (dynamically-consistent)
    # Format: saver.add_dataset("D_adi", (_nsteps, _nadi, _nadi), "C") 
    if "D_adi" in saver.keywords and "D_adi" in saver.np_data.keys():
        saver.save_matrix(t, "D_adi", dm_adi) 

    # Average adiabatic density matrices (raw)
    # Format: saver.add_dataset("D_adi_raw", (_nsteps, _nadi, _nadi), "C") 
    if "D_adi_raw" in saver.keywords and "D_adi_raw" in saver.np_data.keys():
        saver.save_matrix(t, "D_adi_raw", dm_adi_raw) 


    # Average diabatic density matrices (dynamically-consistent)
    # Format: saver.add_dataset("D_dia", (_nsteps, _ndia, _ndia), "C") 
    if "D_dia" in saver.keywords and "D_dia" in saver.np_data.keys():
        saver.save_matrix(t, "D_dia", dm_dia) 

    # Average diabatic density matrices (raw)
    # Format: saver.add_dataset("D_dia_raw", (_nsteps, _ndia, _ndia), "C") 
    if "D_dia_raw" in saver.keywords and "D_dia_raw" in saver.np_data.keys():
        saver.save_matrix(t, "D_dia_raw", dm_dia_raw) 


    # Trajectory-resolved coordinates
    # Format: saver.add_dataset("q", (_nsteps, _ntraj, _dof), "R") 
    if "q" in saver.keywords and "q" in saver.np_data.keys():
        saver.save_matrix(t, "q", q.T()) 

    # Trajectory-resolved momenta
    # Format: saver.add_dataset("p", (_nsteps, _ntraj, _dof), "R") 
    if "p" in saver.keywords and "p" in saver.np_data.keys():
        saver.save_matrix(t, "p", p.T()) 

    # Trajectory-resolved active forces
    # Format: saver.add_dataset("f", (_nsteps, _ntraj, _dof), "R")
#    if "f" in saver.keywords and "f" in saver.np_data.keys():
#        saver.save_matrix(t, "f", f.T())

    # Trajectory-resolved adiabatic TD-SE amplitudes
    # Format: saver.add_dataset("C_adi", (_nsteps, _ntraj, _nadi), "C") 
    if "Cadi" in saver.keywords and "Cadi" in saver.np_data.keys():
        saver.save_matrix(t, "Cadi", Cadi.T()) 

    # Trajectory-resolved diabatic TD-SE amplitudes
    # Format: saver.add_dataset("C_dia", (_nsteps, _ntraj, _ndia), "C") 
    if "Cdia" in saver.keywords and "Cdia" in saver.np_data.keys():
        saver.save_matrix(t, "Cdia", Cdia.T()) 




def save_hdf5_3D_new(saver, i, dyn_var, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type==0:
        t = i


    # Average adiabatic density matrices
    # Format: saver.add_dataset("D_adi", (_nsteps, _nadi, _nadi), "C") 
    if "D_adi" in saver.keywords and "D_adi" in saver.np_data.keys():
        dm_adi = dyn_var.compute_average_dm(1)
        saver.save_matrix(t, "D_adi", dm_adi) 

    # Average diabatic density matrices
    # Format: saver.add_dataset("D_dia", (_nsteps, _ndia, _ndia), "C") 
    if "D_dia" in saver.keywords and "D_dia" in saver.np_data.keys():
        dm_dia = dyn_var.compute_average_dm(0)
        saver.save_matrix(t, "D_dia", dm_dia) 

    # Trajectory-resolved coordinates
    # Format: saver.add_dataset("q", (_nsteps, _ntraj, _dof), "R") 
    if "q" in saver.keywords and "q" in saver.np_data.keys():
        q = dyn_var.get_coords()
        saver.save_matrix(t, "q", q.T()) 

    # Trajectory-resolved momenta
    # Format: saver.add_dataset("p", (_nsteps, _ntraj, _dof), "R") 
    if "p" in saver.keywords and "p" in saver.np_data.keys():
        p = dyn_var.get_momenta()
        saver.save_matrix(t, "p", p.T()) 

    # Trajectory-resolved active forces
    # Format: saver.add_dataset("f", (_nsteps, _ntraj, _dof), "R")
    if "f" in saver.keywords and "f" in saver.np_data.keys():
        f = dyn_var.get_forces()
        saver.save_matrix(t, "f", f.T())

    # Trajectory-resolved adiabatic TD-SE amplitudes
    # Format: saver.add_dataset("C_adi", (_nsteps, _ntraj, _nadi), "C") 
    if "Cadi" in saver.keywords and "Cadi" in saver.np_data.keys():
        Cadi = dyn_var.get_ampl_adi()
        saver.save_matrix(t, "Cadi", Cadi.T()) 

    # Trajectory-resolved diabatic TD-SE amplitudes
    # Format: saver.add_dataset("C_dia", (_nsteps, _ntraj, _ndia), "C") 
    if "Cdia" in saver.keywords and "Cdia" in saver.np_data.keys():
        Cdia = dyn_var.get_ampl_dia()
        saver.save_matrix(t, "Cdia", Cdia.T()) 



def save_hdf5_4D(saver, i, tr, hvib_adi, hvib_dia, St, U, projector, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type==0:
        t = i

    # Trajectory-resolved vibronic Hamiltoninans in the adiabatic representation
    # Format: saver.add_dataset("hvib_adi", (_nsteps, _ntraj, _nadi, _nadi), "C") 
    if "hvib_adi" in saver.keywords and "hvib_adi" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "hvib_adi", hvib_adi) 


    # Trajectory-resolved vibronic Hamiltoninans in the diabatic representation
    # Format: saver.add_dataset("hvib_dia", (_nsteps, _ntraj, _ndia, _ndia), "C") 
    if "hvib_dia" in saver.keywords and "hvib_dia" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "hvib_dia", hvib_dia) 

    # Trajectory-resolved time-overlaps of the adiabatic states
    # Format: saver.add_dataset("St", (_nsteps, _ntraj, _nadi, _nadi), "C") 
    if "St" in saver.keywords and "St" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "St", St) 

    # Trajectory-resolved diabatic-to-adiabatic transformation matrices 
    # Format: saver.add_dataset("basis_transform", (_nsteps, _ntraj, _ndia, _nadi), "C") 
    if "basis_transform" in saver.keywords and "basis_transform" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "basis_transform", U) 

    # Trajectory-resolved projector matrices (from the raw adiabatic to consistent adiabatic)
    # Format: saver.add_dataset("projector", (_nsteps, _ntraj, _nadi, _nadi), "C") 
    if "projector" in saver.keywords and "projector" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "projector", projector) 





def save_tsh_data_123(_savers, params, 
                      i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, states,
                      pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia
                     ):


    hdf5_output_level = params["hdf5_output_level"]
    mem_output_level = params["mem_output_level"]
    txt_output_level = params["txt_output_level"]
    txt2_output_level = params["txt2_output_level"]


    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"]*nsteps)    


    if i%print_freq==0:
        print(F" step= {i}")

    
    if hdf5_output_level>=1 and _savers["hdf5_saver"]!=None:
        save_hdf5_1D(_savers["hdf5_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)

    if mem_output_level>=1 and _savers["mem_saver"]!=None:
        save_hdf5_1D(_savers["mem_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)

    if txt_output_level>=1 and _savers["txt_saver"]!=None:
        save_hdf5_1D(_savers["txt_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)

    if txt2_output_level>=1 and _savers["txt2_saver"]!=None:
        save_hdf5_1D(_savers["txt2_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, 1)



    if hdf5_output_level>=2 and _savers["hdf5_saver"]!=None:
        save_hdf5_2D(_savers["hdf5_saver"], i, states)

    if mem_output_level>=2 and _savers["mem_saver"]!=None:
        save_hdf5_2D(_savers["mem_saver"], i, states)

    if txt_output_level>=2 and _savers["txt_saver"]!=None:
        save_hdf5_2D(_savers["txt_saver"], i, states)

    if txt2_output_level>=2 and _savers["txt2_saver"]!=None:
        save_hdf5_2D(_savers["txt2_saver"], i, states, 1)




    if hdf5_output_level>=3 and _savers["hdf5_saver"]!=None: 
        save_hdf5_3D(_savers["hdf5_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)

    if mem_output_level>=3 and _savers["mem_saver"]!=None: 
        save_hdf5_3D(_savers["mem_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)

    if txt_output_level>=3 and _savers["txt_saver"]!=None: 
        save_hdf5_3D(_savers["txt_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)

    if txt2_output_level>=3 and _savers["txt2_saver"]!=None: 
        save_hdf5_3D(_savers["txt2_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia, 1)





def save_tsh_data_1234_new(_savers, params, i, dyn_var, ham):


    hdf5_output_level = params["hdf5_output_level"]
    mem_output_level = params["mem_output_level"]
    txt_output_level = params["txt_output_level"]
    txt2_output_level = params["txt2_output_level"]

    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"]*nsteps)    


    if i%print_freq==0:
        print(F" step= {i}")

    # =========== Using: def save_hdf5_1D_new(saver, i, params, dyn_var, ham, txt_type=0) =========
    
    if hdf5_output_level>=1 and _savers["hdf5_saver"]!=None:
        save_hdf5_1D_new(_savers["hdf5_saver"], i, params, dyn_var, ham)

    if mem_output_level>=1 and _savers["mem_saver"]!=None:
        #save_hdf5_1D(_savers["mem_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)
        save_hdf5_1D_new(_savers["mem_saver"], i, params, dyn_var, ham)

    if txt_output_level>=1 and _savers["txt_saver"]!=None:
        #save_hdf5_1D(_savers["txt_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)
        save_hdf5_1D_new(_savers["txt_saver"], i, params, dyn_var, ham)

    if txt2_output_level>=1 and _savers["txt2_saver"]!=None:
        #save_hdf5_1D(_savers["txt2_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, 1)
        save_hdf5_1D_new(_savers["txt2_saver"], i, params, dyn_var, ham, 1)


    # =========== Using : def save_hdf5_2D_new(saver, i, dyn_var, txt_type=0) =====================
    if hdf5_output_level>=2 and _savers["hdf5_saver"]!=None:
        #save_hdf5_2D(_savers["hdf5_saver"], i, states)
        save_hdf5_2D_new(_savers["hdf5_saver"], i, dyn_var, ham)

    if mem_output_level>=2 and _savers["mem_saver"]!=None:
        #save_hdf5_2D(_savers["mem_saver"], i, states)
        save_hdf5_2D_new(_savers["mem_saver"], i, dyn_var, ham)

    if txt_output_level>=2 and _savers["txt_saver"]!=None:
        #save_hdf5_2D(_savers["txt_saver"], i, states)
        save_hdf5_2D_new(_savers["txt_saver"], i, dyn_var, ham)

    if txt2_output_level>=2 and _savers["txt2_saver"]!=None:
        #save_hdf5_2D(_savers["txt2_saver"], i, states, 1)
        save_hdf5_2D_new(_savers["txt2_saver"], i, dyn_var, ham, 1)


    # ============ Using: def save_hdf5_3D_new(saver, i, dyn_var, txt_type=0) ==========================

    if hdf5_output_level>=3 and _savers["hdf5_saver"]!=None: 
        #save_hdf5_3D(_savers["hdf5_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)
        save_hdf5_3D_new(_savers["hdf5_saver"], i, dyn_var)

    if mem_output_level>=3 and _savers["mem_saver"]!=None: 
        #save_hdf5_3D(_savers["mem_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)
        save_hdf5_3D_new(_savers["mem_saver"], i, dyn_var)

    if txt_output_level>=3 and _savers["txt_saver"]!=None: 
        #save_hdf5_3D(_savers["txt_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)
        save_hdf5_3D_new(_savers["txt_saver"], i, dyn_var)

    if txt2_output_level>=3 and _savers["txt2_saver"]!=None: 
        #save_hdf5_3D(_savers["txt2_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia, 1)
        save_hdf5_3D_new(_savers["txt2_saver"], i, dyn_var, 1)

    # ============ Using: def save_hdf5_4D_new(saver, i, tr, dyn_var, ham, txt_type=0) ====================

    is_nbra = params["is_nbra"]
    time_overlap_method = params["time_overlap_method"]

    tr_range = []
    if is_nbra ==1:
        tr_range = [0]
    else:
        tr_range = range(dyn_var.ntraj)

    
    for tr in tr_range:

        #if time_overlap_method==0:
        #    x = ham.get_basis_transform(Py2Cpp_int([0, tr]) )            
        #    St = U[tr].H() * x
        #    U[tr] = CMATRIX(x)
        #elif time_overlap_method==1:               
 
        U = ham.get_basis_transform(Py2Cpp_int([0, tr]) )
        St = ham.get_time_overlap_adi(Py2Cpp_int([0, tr]) ) 
        

        if hdf5_output_level>=4 and _savers["hdf5_saver"]!=None: 
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr])) 
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr])) 
            save_hdf5_4D(_savers["hdf5_saver"], i, tr, hvib_adi, hvib_dia, St, U, None)

        if mem_output_level>=4 and _savers["mem_saver"]!=None: 
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr])) 
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr])) 
            save_hdf5_4D(_savers["mem_saver"], i, tr, hvib_adi, hvib_dia, St, U, None)

        if txt_output_level>=4 and _savers["txt_saver"]!=None: 
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr])) 
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr])) 
            save_hdf5_4D(_savers["txt_saver"], i, tr, hvib_adi, hvib_dia, St, U, None)

        if txt2_output_level>=4 and _savers["txt2_saver"]!=None: 
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr])) 
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr])) 
            save_hdf5_4D(_savers["txt2_saver"], i, tr, hvib_adi, hvib_dia, St, U, None, 1)



def print_results12(i, dt, res, prefix, file_output_level):
    """
    Print out the properties at given timestep i to the files

    Args:
        i ( int ): index of the timestep
        dt ( double ): time step size [ units: a.u. ]
        res ( tuple ): contains the properties (at given time) to print out, the order is as follows:

            * q = res[0] ( MATRIX(ndof, ntraj) ): coordiantes of all trajectories
            * p = res[1] ( MATRIX(ndof, ntraj) ): momenta of all trajectories
            * Ekin = res[2] ( double ): average kinetic energy of the ensemble
            * Epot = res[3] ( double ): average potential energy of the ensemble
            * Etot = res[4] ( double ): average total energy of the ensemble
            * dEkin = res[5] ( double ): fluctuation of the average kinetic energy of the ensemble
            * dEpot = res[6] ( double ): fluctuation of the average potential energy of the ensemble
            * dEtot = res[7] ( double ): fluctuation of the average total energy of the ensemble
            * Cadi = res[8] ( CMATRIX(nadi, ntraj) ): amplitudes of adiabatic states for all trajectories
            * Cdia = res[9] ( CMATRIX(ndia, ntraj) ): amplitudes of diabatic states for all trajectories
            * dm_adi = res[10] ( CMATRIX(nadi, nadi) ): adiabatic density matrix averaged over all trajectories
            * dm_dia = res[11] ( CMATRIX(ndia, ndia) ): diabatic density matrix averaged over all trajectories
            * pops = res[12] ( list of nadi doubles ): average populations of all adiabatic states
            * states = res[13] ( list of ntraj ints ): electronic states of all trajectories

        prefix ( string ): the common prefix of the output (usually the directory name)
        file_output_level ( int ): the verbosity of the output

    Returns:
        None

    """

    q, p = res[0], res[1]
    Ekin, Epot, Etot = res[2], res[3], res[4]
    dEkin, dEpot, dEtot = res[5], res[6], res[7]
    Cadi, Cdia = res[8], res[9]
    dm_adi, dm_dia = res[10], res[11]
    pops = res[12]
    states = res[13]
    
    # File output
    if file_output_level >= 1:
        add_doublelist2file("%s/energies.txt" % (prefix), i*dt, [Ekin, Epot, Etot, dEkin, dEpot, dEtot] )
        add_cmatrix2file("%s/D_adi.txt" % (prefix), i*dt, dm_adi )
        add_cmatrix2file("%s/D_dia.txt" % (prefix), i*dt, dm_dia )
        add_matrix2file("%s/SH_pop.txt" % (prefix), i*dt, pops )

    # File output
    if file_output_level >= 2:        
        add_matrix2file("%s/q.txt" % (prefix), i*dt, q )
        add_matrix2file("%s/p.txt" % (prefix), i*dt, p )
        add_cmatrix2file("%s/C_adi.txt" % (prefix), i*dt, Cadi )
        add_cmatrix2file("%s/C_dia.txt" % (prefix), i*dt, Cdia )
        add_intlist2file("%s/states.txt" % (prefix), i*dt, states )        



def print_results3(i, dt, res, prefix, file_output_level, tr):
    """
    Print out the properties at given timestep i to the files

    Args:
        i ( int ): index of the timestep
        dt ( double ): time step size [ units: a.u. ]
        res ( tuple ): contains the properties (at given time) to print out, the order is as follows:

            * hvib_adi = res[0] ( CMATRIX(nadi, nadi) ): adiabatic Hamiltonian for a given trajectory
            * hvib_dia = res[1] ( CMATRIX(ndia, ndia) ): diabatic Hamiltonian for a given trajectory
            * St = res[2] ( CMATRIX(nadi, nadi) ): time-overlap of the adiabatic states for a given trajectory
            * U = res[3] ( CMATRIX(ndia, nadi) ): dia-to-adi transformation matrix for a given trajectory

        prefix ( string ): the common prefix of the output (usually the directory name)
        file_output_level ( int ): the verbosity of the output
        tr ( int ): the index of the trajectory that we handle

    Returns:
        None

    """

    hvib_adi = res[0]
    hvib_dia = res[1]
    St = res[2]
    U = res[3]

    # File output
    if file_output_level >= 3:
        add_cmatrix2file("%s/Hvib_adi_%i.txt" % (prefix, tr), i*dt, hvib_adi )
        add_cmatrix2file("%s/Hvib_dia_%i.txt" % (prefix, tr), i*dt, hvib_dia )
        add_cmatrix2file("%s/St_%i.txt" % (prefix, tr), i*dt, St)
        add_cmatrix2file("%s/basis_transform_%i.txt" % (prefix, tr), i*dt, U )



def read_results(prefix, file_output_level, nadi, ndia, ndof, ntraj):
    """
    This function 
    """
    
    obs_T = None
    obs_q, obs_p = None, None
    obs_Ekin, obs_Epot, obs_Etot = None, None, None
    obs_dEkin, obs_dEpot, obs_dEtot = None, None, None
    obs_Cadi, obs_Cdia = None, None 
    obs_dm_adi, obs_dm_dia = None, None
    obs_pop, obs_states = None, None
    obs_hvib_adi, obs_hvib_dia, obs_St, obs_U = None, None, None, None

    
    if file_output_level >= 1:
        res = data_read.get_data_from_file2("%s/energies.txt" % (prefix), range(0, 7) )
        obs_T = res[0]
        obs_Ekin = res[1]
        obs_Epot = res[2]
        obs_Etot = res[3]
        obs_dEkin = res[4]
        obs_dEpot = res[5]
        obs_dEtot = res[6]
                
        obs_dm_adi = file2cmatrix("%s/D_adi.txt" % (prefix), nadi, nadi)
        obs_dm_dia = file2cmatrix("%s/D_dia.txt" % (prefix), ndia, ndia)
        obs_pop = file2matrix("%s/SH_pop.txt" % (prefix), nadi, 1)
        
    if file_output_level >= 2:
        obs_q = file2matrix("%s/q.txt" % (prefix), ndof, ntraj)
        obs_p = file2matrix("%s/p.txt" % (prefix), ndof, ntraj)
        obs_C_adi = file2cmatrix("%s/C_adi.txt" % (prefix), nadi, ntraj)
        obs_C_dia = file2cmatrix("%s/C_dia.txt" % (prefix), ndia, ntraj)
        obs_states = file2intlist("%s/states.txt" % (prefix))       
        
    if file_output_level >= 3:
        obs_hvib_adi, obs_hvib_dia, obs_St, obs_U = [], [], [], []
        for tr in range(ntraj):                        
            hvib_adi = file2cmatrix("%s/Hvib_adi_%i.txt" % (prefix, tr), nadi, nadi)
            hvib_dia = file2cmatrix("%s/Hvib_dia_%i.txt" % (prefix, tr), nadi, nadi)
            St = file2cmatrix("%s/St_%i.txt" % (prefix, tr), nadi, nadi)
            U = file2cmatrix("%s/basis_transform_%i.txt" % (prefix, tr), ndia, nadi)
            
            obs_hvib_adi.append(hvib_adi)
            obs_hvib_dia.append(hvib_dia)
            obs_St.append(St)
            obs_U.append(U)                


    return obs_T, obs_q, obs_p, obs_Ekin, obs_Epot, obs_Etot, obs_dEkin, obs_dEpot, obs_dEtot, \
           obs_Cadi, obs_Cdia, obs_dm_adi, obs_dm_dia, obs_pop, obs_states, obs_hvib_adi, obs_hvib_dia, obs_St, obs_U

