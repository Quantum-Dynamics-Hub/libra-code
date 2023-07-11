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
       computed/produced during exact on the grid calculations. If you need to print out more 
       data, you shall add the corresponding info in these modules

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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
#import libra_py.units as units
#import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers

#===================== Exact calculations output ====================

def exact_init_hdf5(saver, hdf5_output_level, _nsteps, _ndof, _nstates, _ngrid):

    if hdf5_output_level>=1:

        # Time axis (integer steps)
        saver.add_dataset("timestep", (_nsteps,) , "I")  

        # Time axis
        saver.add_dataset("time", (_nsteps,) , "R")  
        
        # Kinetic energy in diabatic representation
        saver.add_dataset("Ekin_dia", (_nsteps,) , "R")  

        # Kinetic energy in adiabatic representation
        saver.add_dataset("Ekin_adi", (_nsteps,) , "R")  

        # Potential energy in diabatic representation
        saver.add_dataset("Epot_dia", (_nsteps,) , "R")  

        # Potential energy in adiabatic representation
        saver.add_dataset("Epot_adi", (_nsteps,) , "R")  

        # Total energy in diabatic representation
        saver.add_dataset("Etot_dia", (_nsteps,) , "R")  

        # Total energy in adiabatic representation
        saver.add_dataset("Etot_adi", (_nsteps,) , "R")  

        # Wavefunction norm in diabatic representation
        saver.add_dataset("norm_dia", (_nsteps,) , "R")  

        # Wavefunction norm in adiabatic representation
        saver.add_dataset("norm_adi", (_nsteps,) , "R")  


    if hdf5_output_level>=2:

        # Diabatic populations
        saver.add_dataset("pop_dia", (_nsteps, _nstates, 1), "R") 

        # Adiabatic populations
        saver.add_dataset("pop_adi", (_nsteps, _nstates, 1), "R") 


        # All coordinates computed in diabatic rep
        saver.add_dataset("q_dia", (_nsteps, _ndof, 1), "C") 

        # All coordinates computed in adiabatic rep
        saver.add_dataset("q_adi", (_nsteps, _ndof, 1), "C") 

        # All momenta computed in diabatic rep
        saver.add_dataset("p_dia", (_nsteps, _ndof, 1), "C") 

        # All momenta computed in adiabatic rep
        saver.add_dataset("p_adi", (_nsteps, _ndof, 1), "C") 


        # Higher moments computed in diabatic rep
        saver.add_dataset("q2_dia", (_nsteps, _ndof, 1), "C") 

        # Higher moments computed in adiabatic rep
        saver.add_dataset("q2_adi", (_nsteps, _ndof, 1), "C") 

        # Higher moments computed in diabatic rep
        saver.add_dataset("p2_dia", (_nsteps, _ndof, 1), "C") 

        # Higher moments computed in adiabatic rep
        saver.add_dataset("p2_adi", (_nsteps, _ndof, 1), "C") 


    if hdf5_output_level>=3:

        # Density matrix computed in diabatic rep
        saver.add_dataset("denmat_dia", (_nsteps, _nstates, _nstates), "C") 

        # Density matrix computed in adiabatic rep
        saver.add_dataset("denmat_adi", (_nsteps, _nstates, _nstates), "C") 


    if hdf5_output_level>=4:

        # Wavefunction in diabatic rep
        saver.add_dataset("PSI_dia", (_nsteps, _ngrid, _nstates, 1), "C") 

        # Wavefunction in adiabatic rep
        saver.add_dataset("PSI_adi", (_nsteps, _ngrid, _nstates, 1), "C") 

        # Reciprocal wavefunction in diabatic rep
        saver.add_dataset("reciPSI_dia", (_nsteps, _ngrid, _nstates, 1), "C") 

        # Reciprocal wavefunction in adiabatic rep
        saver.add_dataset("reciPSI_adi", (_nsteps, _ngrid, _nstates, 1), "C") 



def exact_init_custom_hdf5(saver, _nsteps, _ncustom_pops, _nstates):

    # Custom populations indx
    saver.add_dataset("custom_pops", (_nsteps, _ncustom_pops, _nstates, 1), "R") 




def save_data_hdf5(step, wfc, saver, params):
    
    dt = params["dt"]    
    masses = Py2Cpp_double(params["masses"])
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level>=1:
        saver.save_scalar(step, "timestep", step) 
        saver.save_scalar(step, "time", step * dt)        
        saver.save_scalar(step, "Ekin_dia", wfc.e_kin(masses, 0)) 
        saver.save_scalar(step, "Ekin_adi", wfc.e_kin(masses, 1))         
        saver.save_scalar(step, "Epot_dia", wfc.e_pot(0)) 
        saver.save_scalar(step, "Epot_adi", wfc.e_pot(1))         
        saver.save_scalar(step, "Etot_dia", wfc.e_tot(masses, 0)) 
        saver.save_scalar(step, "Etot_adi", wfc.e_tot(masses, 1))                 
        saver.save_scalar(step, "norm_dia", wfc.norm(0) )
        saver.save_scalar(step, "norm_adi", wfc.norm(1) )
        
        
    if hdf5_output_level>=2:    

        saver.save_matrix(step, "pop_dia", wfc.get_pops(0) ) 
        saver.save_matrix(step, "pop_adi", wfc.get_pops(1) )     
        saver.save_matrix(step, "q_dia", wfc.get_pow_q(0, 1) ) 
        saver.save_matrix(step, "q_adi", wfc.get_pow_q(1, 1) )         
        saver.save_matrix(step, "q2_dia", wfc.get_pow_q(0, 2) ) 
        saver.save_matrix(step, "q2_adi", wfc.get_pow_q(1, 2) )         
        saver.save_matrix(step, "p_dia", wfc.get_pow_p(0, 1) ) 
        saver.save_matrix(step, "p_adi", wfc.get_pow_p(1, 1) )         
        saver.save_matrix(step, "p2_dia", wfc.get_pow_p(0, 2) ) 
        saver.save_matrix(step, "p2_adi", wfc.get_pow_p(1, 2) ) 

        if "custom_pops" in params.keys():
            ncustom_pops = len(params["custom_pops"])
            for pop_type in range(ncustom_pops):

                pop_prms = params["custom_pops"][pop_type]            
                pop_val = wfc.get_pops(pop_prms[0], Py2Cpp_double(pop_prms[1]), Py2Cpp_double(pop_prms[2])) 

                saver.save_multi_matrix(step, pop_type, "custom_pops",  pop_val )

        
        
    if hdf5_output_level>=3:                
        saver.save_matrix(step, "denmat_dia", wfc.get_den_mat(0) ) 
        saver.save_matrix(step, "denmat_adi", wfc.get_den_mat(1) ) 
        
        
    if hdf5_output_level>=4:                
        """
        This is a REALLY-REALLY bad option to go since it is soo time-consuming
        """
        
        for ipt in range(wfc.Npts):
            saver.save_multi_matrix(step, ipt, "PSI_dia", wfc.PSI_dia[ipt]) 
            
        for ipt in range(wfc.Npts):            
            saver.save_multi_matrix(step, ipt, "PSI_adi", wfc.PSI_adi[ipt]) 
            
        for ipt in range(wfc.Npts):
            saver.save_multi_matrix(step, ipt, "reciPSI_dia", wfc.reciPSI_dia[ipt]) 
            
        for ipt in range(wfc.Npts):            
            saver.save_multi_matrix(step, ipt, "reciPSI_adi", wfc.reciPSI_adi[ipt])             


            

def save_data_mem(step, wfc, saver, params):
    
    dt = params["dt"]    
    masses = Py2Cpp_double(params["masses"])
    mem_output_level = params["mem_output_level"]
    
    
    if mem_output_level>=1:         
        saver.add_data("timestep", step)
        saver.add_data("time", step * dt)
        saver.add_data("Ekin_dia", wfc.e_kin(masses, 0) )
        saver.add_data("Ekin_adi", wfc.e_kin(masses, 1) )
        saver.add_data("Epot_dia", wfc.e_pot(0) )
        saver.add_data("Epot_adi", wfc.e_pot(1) )
        saver.add_data("Etot_dia", wfc.e_tot(masses, 0) )
        saver.add_data("Etot_adi", wfc.e_tot(masses, 1) )
        saver.add_data("norm_dia", wfc.norm(0) )
        saver.add_data("norm_adi", wfc.norm(1) )

    if mem_output_level>=2:        
        saver.add_data("pop_dia", wfc.get_pops(0) ) 
        saver.add_data("pop_adi", wfc.get_pops(1) )     
        saver.add_data("q_dia", wfc.get_pow_q(0, 1) ) 
        saver.add_data("q_adi", wfc.get_pow_q(1, 1) )         
        saver.add_data("q2_dia", wfc.get_pow_q(0, 2) ) 
        saver.add_data("q2_adi", wfc.get_pow_q(1, 2) )         
        saver.add_data("p_dia", wfc.get_pow_p(0, 1) ) 
        saver.add_data("p_adi", wfc.get_pow_p(1, 1) )         
        saver.add_data("p2_dia", wfc.get_pow_p(0, 2) ) 
        saver.add_data("p2_adi", wfc.get_pow_p(1, 2) )         

    if mem_output_level>=3:        
        saver.add_data("denmat_dia", wfc.get_den_mat(0) )  
        saver.add_data("denmat_adi", wfc.get_den_mat(1) )         
        
    if mem_output_level>=4:                
        saver.add_data("PSI_dia", wfc.PSI_dia) 
        saver.add_data("PSI_adi", wfc.PSI_adi) 
        saver.add_data("reciPSI_dia", wfc.reciPSI_dia) 
        saver.add_data("reciPSI_adi", wfc.reciPSI_adi) 
        


def init_tsh_savers(params, model_params, nsteps, wfc):


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
        exact_init_hdf5(_savers["hdf5_saver"], hdf5_output_level, nsteps, wfc.ndof, wfc.nstates, wfc.Npts)


    #====== TXT ========
    txt_output_level = params["txt_output_level"]
    if params["txt_output_level"] > 0:
        _savers["txt_saver"] = data_savers.mem_saver(properties_to_save)
        #_savers["txt_saver"].set_compression_level(params["use_compression"], params["compression_level"])


    #====== TXT2: No intermediate memory allocation ========
    txt2_output_level = params["txt2_output_level"]

    if params["txt2_output_level"] > 0:
        _savers["txt2_saver"] = data_savers.mem_saver(properties_to_save)

        # Here, nsteps is set to 1 since in this type of saver, we only care about the current values,
        # not all the timesteps - that would be too consuming
        #exact_init_hdf5(_saver["hdf5_save"], txt2_output_level, 1, wfc.ndof, wfc.nstates, wfc.Npts)


    #====== MEM =========
    mem_output_level = params["mem_output_level"]

    if mem_output_level > 0:
        _savers["mem_saver"] =  data_savers.mem_saver(properties_to_save)
        #print("Saver not implemented")
        #sys.exit(0)
        #init_tsh_data(_savers["mem_saver"], mem_output_level, nsteps, ntraj, nnucl, nadi, ndia)


    return _savers



