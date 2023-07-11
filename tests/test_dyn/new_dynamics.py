#*********************************************************************************
#* Copyright (C) 2022 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""

 This is a script to run the NAMD with model Hamiltonians using the newest version
 of the Libra's dynamics

"""

import sys
import cmath
import math
import os
import h5py
import matplotlib.pyplot as plt   # plots
import numpy as np
import time
import warnings

from liblibra_core import *
import util.libutil as comn
from libra_py import units
import libra_py.models.Holstein as Holstein
import libra_py.models.Morse as Morse
from libra_py import dynamics_plotting
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot



#from matplotlib.mlab import griddata
#%matplotlib inline 
warnings.filterwarnings('ignore')

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





def compute_model(q, params, full_id):

    model = params["model"]
    res = None
    
    if model==1:        
        res = Holstein.Holstein2(q, params, full_id) 
    elif model==2:
        pass #res = compute_model_nbra(q, params, full_id)
    elif model==3:
        pass #res = compute_model_nbra_files(q, params, full_id)        
    elif model==4:
        res = Morse.general(q, params, full_id)    

    return res




model_params1 = {"model":1, "model0":1, "E_n":[0.0,  0.0], "x_n":[0.0,  2.5],"k_n":[0.002, 0.005],"V":0.000, "nstates":2}
model_params2 = {"model":1, "model0":1, "E_n":[0.0,  0.0], "x_n":[0.0,  2.5],"k_n":[0.002, 0.005],"V":0.001, "nstates":2}
model_params3 = {"model":1, "model0":1, "E_n":[0.0,  0.0], "x_n":[0.0,  2.5],"k_n":[0.002, 0.005],"V":0.01, "nstates":2}



decoherence_rates = MATRIX(2,2)
ave_gaps = MATRIX(2,2)
schwartz_decoherence_inv_alpha = MATRIX(2,2)


dyn_general = { "nsteps":2500, "dt":10.0, "ntraj":1, 
                "nac_update_method":2, "hvib_update_method":1,
                "force_method":1, "rep_force":1,
                "use_boltz_factor":0, "num_electronic_substeps":1,
                "nstates":2, "which_adi_states":range(2), "which_dia_states":range(2),
                "state_tracking_algo":2, "do_phase_correction":1,
                "num_electronic_substeps":1, "isNBRA":0, "is_nbra":0, "progress_frequency":1.0,
                "decoherence_rates":decoherence_rates, "ave_gaps":ave_gaps,
                "schwartz_decoherence_inv_alpha":schwartz_decoherence_inv_alpha,

                "tsh_method":-1, "mem_output_level":3,
                "properties_to_save":[ "timestep", "time", "q", "p", "Cadi", "Cdia", "Epot_ave", "Ekin_ave", "Etot_ave",
                "se_pop_adi", "se_pop_dia", "sh_pop_adi" ],
                "prefix":"adiabatic_md", "prefix2":"adiabatic_md"
              }
rnd = Random()

model_params = model_params1



#=================== Dynamics =======================

#elec_params = {"ndia":2, "nadi":2, "init_type":3, "rep":1, "istates":[1.0, 0.0], "verbosity":1, "init_dm_type":0 }
nucl_params = {"ndof":1, "init_type":2, "q":[-4.0], "p":[0.0], "mass":[2000.0], "force_constant":[0.01] }
elec_params = {}
dyn_params = dict(dyn_general)

case = 3

if case==0:  # adiabatic init, adiabatic TDSE
    elec_params = {"ndia":2, "nadi":2, "init_type":3, "rep":1, "istates":[1.0, 0.0], "verbosity":1, "init_dm_type":0 , "verbosity":2 }
    dyn_params.update({"rep_tdse":1, "ham_update_method":1, "ham_transform_method":1, "tsh_method":0, "ntraj":25, "nsteps":2500, 
                       "prefix":F"case{case}", "prefix2":F"case{case}" })  

elif case==1: # diabatic init, adiabatic TDSE
    elec_params = {"ndia":2, "nadi":2, "init_type":3, "rep":0, "istates":[1.0, 0.0], "verbosity":1, "init_dm_type":0 , "verbosity":2 }
    dyn_params.update({"rep_tdse":1, "ham_update_method":1, "ham_transform_method":1, "tsh_method":0, "ntraj":25, "nsteps":2500,
                        "prefix":F"case{case}", "prefix2":F"case{case}" })  

elif case==2:  # adiabatic init, diabatic TDSE
    elec_params = {"ndia":2, "nadi":2, "init_type":3, "rep":1, "istates":[1.0, 0.0], "verbosity":1, "init_dm_type":0 , "verbosity":2 }
    dyn_params.update({"rep_tdse":0, "ham_update_method":1, "ham_transform_method":1, "tsh_method":-1, "ntraj":25, "nsteps":2500,
                       "prefix":F"case{case}", "prefix2":F"case{case}" })  

elif case==3:  # diabatic init, diabatic TDSE
    elec_params = {"ndia":2, "nadi":2, "init_type":3, "rep":0, "istates":[1.0, 0.0], "verbosity":1, "init_dm_type":0 , "verbosity":2 }
    dyn_params.update({"rep_tdse":0, "ham_update_method":1, "ham_transform_method":1, "tsh_method":-1, "ntraj":25, "nsteps":2500,
                       "prefix":F"case{case}", "prefix2":F"case{case}" })  



res = tsh_dynamics.generic_recipe(dyn_params, compute_model, model_params, elec_params, nucl_params, rnd)


#================== Surfaces =================
plot_params = {"colors": colors, "clrs_index": clrs_index, "ylim":[-0.01, 0.06], "xlim":[-4, 5], "prefix":F"case{case}", "save_figures":1, "do_show":0 }
dynamics_plotting.plot_surfaces(compute_model, [ model_params ], [0, 1], -4.0, 5.0, 0.05, plot_params)


#============ Plotting ==================


plot_params = { "prefix":F"case{case}", "filename":"mem_data.hdf", "output_level":3,
                "which_trajectories":[0], "which_dofs":[0], "which_adi_states":[0,1], "which_dia_states":[0,1], 
                "frameon":True, "linewidth":5, "dpi":300,
                "axes_label_fontsize":(8,8), "legend_fontsize":8, "axes_fontsize":(8,8), "title_fontsize":8,
                "what_to_plot":["coordinates", "momenta", "energies", "phase_space", "se_pop_adi", "se_pop_dia", "sh_pop_adi"], 
                "which_energies":["potential", "kinetic", "total"],
                "save_figures":1, "do_show":0
              }

tsh_dynamics_plot.plot_dynamics(plot_params)
