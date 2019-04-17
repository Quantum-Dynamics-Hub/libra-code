#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file defaults.py
# This module implements the functions that set up default parameters for different
# types of simulations

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


def set_defaults(params, interface, recipe=""):
## 
# \param[in, out] params (dictionary) A dictionary of the parameters for simulaions. You pass an empty
# (or not) dictionary to this function, and it will populate the dictionary with the default
# values, depending on type of interface you about to use (controlled by the "interface" argument)
# and the type of simulation you want to use (controlled by the "recipe")
# \param[in] interface (string) A type of Libra-X interface you are going to use in the calculations
# Options are: "QE", "GAMESS"
# \param[in] recipe (string) The flag controlling specific selection of settings for a specific
# type of simulations (recipe) you want to perform.
# Options (yet to be defined): ""

    if interface=="GAMESS" or interface=="QE" or interface=="G09":
        params["interface"] = interface
        print "Setting default parameters for Libra-",interface, " calculations"
    else:
        print "Error in set_defaults: The interface ", interface, " is not known"
        sys.exit(0)

    # Set up general parameters for both QE and GAMESS runs


    ##### select directories where the results will be printed out. #####
    cwd = os.getcwd()

    # TO DO: for the following 3 lines of code, we'll need to verify if the directories
    # exist, and create them if they are not    
    # directory where the all results will be printed out
    params["res"] = cwd+"/res/" 

    # directory where MO basis vibronic Hamiltonians will be printed out
    params["mo_ham"] = cwd+"/mo_ham/" 

    # directory where SD basis vibronic hamiltonians will be printed out
    params["sd_ham"] = cwd+"/sd_ham/" 

    # directory where hop probabilities will be printed out
    params["hop_probs"] = cwd+"/hop_probs/"

    if interface=="QE":
        # The file into which the MD trajectorie will be printed out
        params["traj_file_prefix"] = params["res"]+"md"  

        # The file into which the kinetic, potential, system, and 
        # the extended (system + thermostat) energies will be printed out
        params["ene_file_prefix"] = params["res"]+"ene"

        # The file into which dipole moment matrices will be printed out
        params["mu_file_prefix"] = params["res"]+"mu"     
 
        # The prefix for the files containing the SE population 
        # (if velocity rescaling is applied, this is averaged over stochastic realizations
        # of the TSH processes). 
        # The file names are defined as se_pop_"+initial geometry+"_"initial excitation",
        # where initial_geometry is an index of the given starting geometry (positions,velocities)
        params["se_pop_file_prefix"] = params["res"]+"se_pop"

        # Similar to the above, but now for SH-derived populations
        params["sh_pop_file_prefix"] = params["res"]+"sh_pop"

        # The prefix for the files containing the SE population averaged over initial geometries.
        # The file names are derives as "se_pop_ex"+initial excitation, where
        # initial_excitation is an index of particular starting excitaion
        params["se_pop_ex_file_prefix"] = params["res"]+"se_pop_ex" 

        # Similar to the above, but for the SH-derived populations
        params["sh_pop_ex_file_prefix"] = params["res"]+"sh_pop_ex"


    ##### Flags for debugging and print out control #####
    # Whether to print auxiliary results. A large amount of printout (MD, energy 
    # trajectories, etc..) will be produced, so use this only when necessary. 
    # For now, it is set to 1 since we still develop the code
    # Options: 0 -> no, 1 -> yes
    params["print_aux_results"] = 1 

    # Whether to compute and print electronic coherences (c^*_i * c_j)
    # Options: 0 -> no, 1 -> yes
    params["print_coherences"] = 0     

    # Print SD basis vibronic Hamiltonian
    # Options: 0 -> no, 1 -> yes
    params["print_sd_ham"] = 0 

    # Print full and reduced size MO basis vibronic Hamiltonian
    # Options: 0 -> no, 1 -> yes
    params["print_mo_ham"] = 0

    # Print the debug info into standard output: hopping probabilities matrices and SH_states 
    # Options: 0 -> no, 1 -> yes
    params["print_tsh_probabilities"] = 1 

    # Print the hopping probabilities if they are larger than 1.
    # (To check whether dt_nucl is too large or not.)
    params["check_tsh_probabilities"] = 1 

    # Flag for printing out the information about atomic orbital basis
    # Options: 1 -> yes. otherwise -> no. 
    # Don't choose 1 when you use PM6: PM6 calculation doesn't output it at present.
    if interface=="GAMESS":
        params["flag_ao"] = 1
    elif interface=="QE":
        params["flag_ao"] = 0


    if interface=="GAMESS":
        # Flag to print MD, Energy, and dipole moment results of 
        # the SH calculation with velocity rescaling
        # Options: 0 -> no, 1 -> yes
        params["print_SH_results_with_scaling"] = 0 

        # Flag to print the debug info into standard output: density matrices
        # for the wavefunctions at different time steps
        # Options: 0 -> no, 1 -> yes
        params["debug_densmat_output"] = 1  

        # Flag to print the debug info into standard output: transition dipole 
        # moment matrices
        # Options: 0 -> no, 1 -> yes
        params["debug_mu_output"] = 0 

        # Flag to print the debug info into standard output: unpacked data
        # Options: 0 -> no, 1 -> yes
        params["debug_gms_unpack"] = 0  

    elif interface=="QE":
        params["qe_debug_print"] = 0


    ##### General simulation (methodology definition) parameters #####
    
    # An algorithm for surface hopping 
    # Options:  1 -> FSSH, 2 -> GFSH , 3 -> MSSH
    params["tsh_method"] = 1  

    # Representation - this parameter is used only to determine
    # how to perform velocity rescaling
    # Options: 0 -> diabatic (uniform rescaling, needs only energies, not derivative couplings)
    #          1 - adiabatic - derivative coupling vectors are needed
    params["rep"] = 0      

    # A flag to select the Boltzmann scaling in lieu of the hop rejection and 
    # velocity rescaling
    # Options: 0 -> do not use Boltzmann scaling (in case you do velocity rescaling),
    #          1 -> use Boltzmann scaling (CPA)
    params["use_boltz_factor"] = 0    

    # The flag to control velocity rescaling.
    # Options: 0 -> no velocity rescaling, 1 -> do rescaling
    params["do_rescaling"] = 1           

    # The option that determines what to do if the hop was rejected because
    #  of the energy conservation(frustrated hop)
    # Options: 0 -> do not change the velocity, 1 - reverse its direction
    params["do_reverse"] = 1

    # The option that depends on how we assume the basis states are
    # Options: 0 -> orthogonal (good for GAMESS and Pyxaid-type wavefunction)
    #          1 -> non-orthogonal (as in QE with delta-SCF)
    params["non-orth"] = 0

    # The option that determines if decoherence effects are included after hops. 
    # Options: 0 -> no decoherence
    #          1 -> decoherence
    params["do_collapse"] = 0

    # Whether to use any thermostat
    # Options: None - not to use
    params["therm"] = None

    ##### Optional MM interactions on top of the QM #####

    # The flag controlling whether to include MM interactions
    # Options: 0 -> no, 1 -> yes
    params["is_MM"] = 0 

    # The fraction of the MM (QM) part in the QM/MM mixing:
    # E_total = qm_frac*E(QM) + mm_frac*E(MM), same for forces (at least for the ground state)!
    # Consider: qm_frac = 1, mm_frac = 1 and ff with only vdw interactions - that would be
    #                        an empirically-corrected dispersion for the QM part
    #  or       qm_frac = 0, mm_frac = 1 and ff with all types of interactions - that would be
    #                        just a classical dynamics (note: the QM calculations may still be done
    #                        e.g. to compute couplings etc.
    params["MM_fraction"] = 0.0
    params["QM_fraction"] = 1.0

    # create a Universe object
    params["U"] = Universe()
    LoadPT.Load_PT(params["U"], "elements.txt")

    # Create a Force field object
    params["ff"] = ForceField({"mb_functional":"LJ_Coulomb","R_vdw_on": 10.0,"R_vdw_off":15.0 })
    LoadUFF.Load_UFF(params["ff"], "uff.d")

    # The file contaiining the atomic connectivity information for the MM part
    params["ent_file"] = ""
 


