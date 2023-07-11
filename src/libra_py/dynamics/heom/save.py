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
       computed/produced during HEOM calculations. If you need to print out more 
       data, you shall add the corresponding info in these modules

       List of functions:
           * init_heom_data(saver, hdf5_output_level, _nsteps, _nquant)
           * init_heom_savers(params, nquant)
           * save_heom_hdf5(step, saver, params, denmat)
           * save_heom_data(_savers, step, print_freq, params, rho_unpacked)

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






def init_heom_savers(params, nquant):

    #================ Create savers ==================    
    prefix = params["prefix"]

    # Create an output directory, if not present    
    if not os.path.isdir(prefix):
        os.mkdir(prefix)

    properties_to_save = params["properties_to_save"]


    _savers = {"hdf5_saver":None, "txt_saver":None, "mem_saver":None }

    #====== HDF5 ========
    hdf5_output_level = params["hdf5_output_level"]
    
    if hdf5_output_level > 0:                
        _savers["hdf5_saver"] = data_savers.hdf5_saver(F"{prefix}/data.hdf", properties_to_save) 
        _savers["hdf5_saver"].set_compression_level(params["use_compression"], params["compression_level"])
        init_heom_data(_savers["hdf5_saver"], hdf5_output_level, params["nsteps"], nquant)

    #====== TXT ========
    if params["txt_output_level"] > 0:
        pass
    
    #====== MEM =========
    mem_output_level = params["mem_output_level"]

    if mem_output_level > 0:
        _savers["mem_saver"] =  data_savers.mem_saver(properties_to_save)
        init_heom_data(_savers["mem_saver"], mem_output_level, params["nsteps"], nquant)

    return _savers                         
    





def save_heom_hdf5(step, saver, params, denmat):
    
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



def save_heom_data(_savers, step, print_freq, params, rho_unpacked):

    #================ Saving the data ==================

    #if step%print_freq==0:
    #    print(F" step= {step}")
        
    # Save properties
    if _savers["hdf5_saver"] != None:            
        save_heom_hdf5(step, _savers["hdf5_saver"], params, rho_unpacked[0])
        
    if _savers["txt_saver"] != None:            
        pass
        
    if _savers["mem_saver"] != None:            
        prms = dict(params)
        prms["hdf5_output_level"] = prms["mem_output_level"]
        save_heom_hdf5(step, _savers["mem_saver"], prms, rho_unpacked[0])

    
