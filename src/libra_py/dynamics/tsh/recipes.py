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
.. module:: recipes
   :platform: Unix, Windows
   :synopsis: This module implements a number of predefined simulation parameters sets, for convenience

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

from libra_py import units

def adiabatic_md_interfaces_params():
    """
    This function returns the generic parameters for adiabatic MD. 
    It assumes everything is done in the adiabatic basis, and the minimal number of adiabatic 
    properties is available (energies, gradients). It is most suitable for the runs that rely on
    interfaces of Libra with external electronic structure packages. No decoherence, no phase correction,
    no NACs, no surface hopping, no trajectory entanglement.
    """


    # Dynamics (Simulation) parameters
    dyn_params = {}

    # 1 - solve the TDSE for electronic amplitudes in adiabatic representation
    dyn_params["rep_tdse"] = 1

    # 1 - Hamiltonian is computed directly in the adiabatic representation
    dyn_params["rep_ham"] = 1

    # 1 - use adiabatic properties for SH (if any)
    dyn_params["rep_sh"] = 1

    # 1 - use adiabatic properties for LZ SH (not relevant for this case)
    dyn_params["rep_lz"] = 1

    # -1 - adibatic MD, consider:  0 - FSSH, 1- GFSH, 2- MSSH, 3 - DISH
    dyn_params["tsh_method"] = -1

    # 1 - use state-specific forces
    dyn_params["force_method"] = 1

    # 0 - don't update NACs and vibronic Ham: we'd need derivative couplings
    dyn_params["nac_update_method"] = 0

    # 1 - forces are computed in the adiabatic representation
    dyn_params["rep_force"] = 1

    # Options to control the acceptance of the proposed hops:
    # 0 - accept all; consider: 
    # 10 - adiabatic energies
    # 21 - difference of state-specific forces
    # 31 - quantum Boltzmann factors (e.g. NBRA)
    dyn_params["hop_acceptance_algo"] = 0 

    # Options to control momenta changes upon successful or frustrated hops:
    # 0 - don't rescale; consider: 
    # 100 - based on adiabatic energy, don't reverse on frustrated hops                          
    # 101 - based on adiabatic energy, reverse on frustrated hops
    dyn_params["momenta_rescaling_algo"] = 0      
  
    # Whether to scale the SH probabilities by the Boltzmann factor: 0 - do not scale, 1 - scale
    dyn_params["use_boltz_factor"] = 0      

    # Temperature, K
    dyn_params["Temperature"] = 300.0

    # Option to perform the phase correction: 0 - no, 1 - yes (default) 
    dyn_params["do_phase_correciton"] = 0

    # the minimal value of the time-overlap for considering 
    # phase corrections - no correction applied if the time-overlap is smaller in magnitude than 
    # this parameter
    dyn_params["tol"] = 1e-3 

    # state tracking algo: 0 - no, 1 - Sato, 2 - mincost (Munkres-Kuhn, MK)
    dyn_params["state_tracking_algo"] = 0

    # parameters to control the mincost state trancking algo (irrelevant here)
    dyn_params["MK_alpha"] = 0.0
    dyn_params["MK_verbosity"] = 0

    # 0 - don't couple (ntangle) trajectories
    dyn_params["entanglement_opt"] = 0

    # ETHD3 (entangled trajectories) options - irrelevant here
    dyn_params["ETHD3_alpha"] = 0.0
    dyn_params["ETHD3_beta"] = 0.0

    # decoherence algo: -1 - no, 0 - MSDM and alike, 1 - ID-A/ID-S/ID-C
    dyn_params["decoherence_algo"] = -1

    # Matrix of the decoherence rates (used with some methods) - irrelevant here
    #dyn_params["decoherence_rates"] = DR

    # Type of dephasing times/rates calculation (only if decoherence is used)
    # 0 - use the rates read out from the input  [default]
    # 1 - use the energy-based decoherence method (EDC)   
    dyn_params["decoherence_times_type"] = 0 

    # An empirical parameter used in the EDC method: [default = 1.0 Ha]
    dyn_params["decoherence_C_param"] = 1.0 

    # An empirical parameter used in the EDC method: [default = 0.1 Ha]
    dyn_params["decoherence_eps_param"] = 1.0 

    # A flag to apply the dephasing-informed approach of Sifain et al 
    # to correct dephasing times: 0 - don't apply; 1 - use it 
    dyn_params["dephasing_informed"] = 0
  
    # A matrix that contains the averaged moduli of the energy gaps:
    # E_ij = <|E_i - E_j|>
    # It is needed when dephasing_informed option is used
    #dyn_params["ave_gaps"] = AG
  
    # Option to control the instantaneous decoherence methodology,
    # only used with decoherence_algo == 1
    # 0 - ID-S; 1 - ID-A [default]; 2 - ID-C - consistent ID - an experimental algorithm
    dyn_params["instantaneous_decoherence_variant"] = 1

    # How to collapse wavefunction amplitudes in the decoherence schemes:
    # 0 - by rescaling the magnitude of the amplitude vector elements, but preserving "phase" [ default ]
    # 1 - by resetting the amplitudes to 1.0+0.0j. This option changes phase 
    dyn_params["collapse_option"] = 0
  
    # integration timestep [units: a.u., default: 41 a.u. = 1 fs]
    dyn_params["dt"] = 1.0 * units.fs2au

    # how many nuclear dynamics steps to perform
    dyn_params["nsteps"] = 1

    # Ensemble: which ensemble to use: 0 - NVE, 1 - NVT
    dyn_params["ensemble"] = 0 

    # thermostat, only if NVT is used
    dyn_params["thermostat_params"] = { }
  

    # What to save
    dyn_params["properties_to_save"] = [ "timestep", "time", "Ekin_ave", "Epot_ave", "Etot_ave",      
                                         "dEkin_ave", "dEpot_ave", "dEtot_ave", "states", "SH_pop", "SH_pop_raw",
                                         "q", "p", "Cadi", "Cdia", "hvib_adi",
                                         "D_adi", "D_adi_raw", "D_dia", "D_dia_raw", 
                                         "hvib_adi", "hvib_dia",  "St", "basis_transform", "projector"   ] 
    # How to save
    dyn_params.update( { "hdf5_output_level":-1, "prefix":"out", 
                         "mem_output_level":3,
                         "use_compression":0, "compression_level":[0,0,0], 
                         "txt_output_level":0,
                        "progress_frequency":1.0
                       } 
                      )

    return dyn_params
