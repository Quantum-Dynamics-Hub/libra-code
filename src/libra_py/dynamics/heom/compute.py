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
.. module:: compute
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
import time

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


import util.libutil as comn
import libra_py.units as units
from . import save


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
                       "num_threads":1
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




    # Initialize savers
    _savers = save.init_heom_savers(params, nquant) 

    #============== Propagation =============

    unpack_rho(rho_unpacked, rho)
    start = time.time()
    for step in range(params["nsteps"]):

        save.save_heom_data(_savers, step, print_freq, params, rho_unpacked)


        if step % params["filter_after_steps"] == 1:
            unpack_rho(rho_unpacked, rho)
            params["zero"] = filter(rho_unpacked, params["tolerance"]);
            pack_rho(rho_unpacked, rho)

        rho = RK4(rho, params["dt"], compute_heom_derivatives, params)

    end = time.time()
    print(F"Calculations took {end - start} seconds")


    # For the mem_saver - store all the results into HDF5 format only at the end of the simulation    
    if _savers["mem_saver"] != None:        
        prefix = params["prefix"]
        _savers["mem_saver"].save_data( F"{prefix}/mem_data.hdf", params["properties_to_save"], "w")
        return _savers["mem_saver"]




