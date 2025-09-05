# *********************************************************************************
# * Copyright (C) 2019-2020 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
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

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

import util.libutil as comn
# import libra_py.units as units
# import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers


# ===================== TSH calculations output ====================

def init_tsh_data(saver, output_level, _nsteps, _ntraj, _ndof, _nadi, _ndia):
    """
    saver - can be either hdf5_saver or mem_saver

    """
    # print(F"In init_tsh_data with:nsteps = {_nsteps}, ntraj = {_ntraj}, ndof = {_ndof}, nadi = {_nadi}, ndia = {_ndia} ")
    # print(F"output_level = {output_level}")

    if output_level >= 1:

        # Time axis (integer steps)
        if "timestep" in saver.keywords:
            saver.add_dataset("timestep", (_nsteps,), "I")

        # Time axis
        if "time" in saver.keywords:
            saver.add_dataset("time", (_nsteps,), "R")

        # Average kinetic energy
        if "Ekin_ave" in saver.keywords:
            saver.add_dataset("Ekin_ave", (_nsteps,), "R")

        # Average potential energy
        if "Epot_ave" in saver.keywords:
            saver.add_dataset("Epot_ave", (_nsteps,), "R")

        # Average total energy
        if "Etot_ave" in saver.keywords:
            saver.add_dataset("Etot_ave", (_nsteps,), "R")

        # Fluctuation of average kinetic energy
        if "dEkin_ave" in saver.keywords:
            saver.add_dataset("dEkin_ave", (_nsteps,), "R")

        # Fluctuation of average potential energy
        if "dEpot_ave" in saver.keywords:
            saver.add_dataset("dEpot_ave", (_nsteps,), "R")

        # Fluctuation of average total energy
        if "dEtot_ave" in saver.keywords:
            saver.add_dataset("dEtot_ave", (_nsteps,), "R")

        # Thermostat energy
        if "Etherm" in saver.keywords:
            saver.add_dataset("Etherm", (_nsteps,), "R")

        # System + thermostat energy
        if "E_NHC" in saver.keywords:
            saver.add_dataset("E_NHC", (_nsteps,), "R")

        # TC-NBRA kinetic energy - averaged over all trajectories
        if "tcnbra_ekin" in saver.keywords:
            saver.add_dataset("tcnbra_ekin", (_nsteps,), "R")

        # TC-NBRA thermostat energy
        if "tcnbra_thermostat_energy" in saver.keywords:
            saver.add_dataset("tcnbra_thermostat_energy", (_nsteps,), "R")

        # Average kinetic energy in QTSH
        if "Ekin_ave_qtsh" in saver.keywords:
            saver.add_dataset("Ekin_ave_qtsh", (_nsteps,), "R")

        # KC-RPMD auxiliary electronic variable kinetic energy
        if "ekin_aux_var" in saver.keywords:
            saver.add_dataset("ekin_aux_var", (_nsteps,), "R")

    if output_level >= 2:

        # Trajectory-resolved instantaneous adiabatic states
        if "states" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("states", (_nsteps, _ntraj), "I")

        # Trajectory-resolved instantaneous adiabatic states
        if "states_dia" in saver.keywords:  # and "states_dia" in saver.np_data.keys():
            saver.add_dataset("states_dia", (_nsteps, _ntraj), "I")

        # Average adiabatic SE populations
        if "se_pop_adi" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("se_pop_adi", (_nsteps, _nadi), "R")

        # Average diabatic SE populations
        if "se_pop_dia" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("se_pop_dia", (_nsteps, _ndia), "R")

        # Average adiabatic SH populations
        if "sh_pop_adi" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("sh_pop_adi", (_nsteps, _nadi), "R")

        # Average diabatic SH populations
        if "sh_pop_dia" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("sh_pop_dia", (_nsteps, _ndia), "R")

        # Average adiabatic SH populations from the Tempelaar and Reichman's method
        if "sh_pop_adi_TR" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("sh_pop_adi_TR", (_nsteps, _nadi), "R")

        # Average diabatic SH populations from the Tempelaar and Reichman's method
        if "sh_pop_dia_TR" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("sh_pop_dia_TR", (_nsteps, _nadi), "R")

        # Average adiabatic MASH populations
        if "mash_pop_adi" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("mash_pop_adi", (_nsteps, _nadi), "R")

        # Average diabatic MASH populations
        if "mash_pop_dia" in saver.keywords:  # and "states" in saver.np_data.keys():
            saver.add_dataset("mash_pop_dia", (_nsteps, _ndia), "R")

        # Average errors from FSSH3
        if "fssh3_average_errors" in saver.keywords:  # and "fssh3_average_errors" in saver.np_data.keys():
            saver.add_dataset("fssh3_average_errors", (_nsteps, 5), "R")

        # Trajectory-resolved auxiliary electronic variable coordinate
        if "y_aux_var" in saver.keywords:  # and "y_aux_var" in saver.np_data.keys():
            saver.add_dataset("y_aux_var", (_nsteps, 1), "R")

        # Trajectory-resolved auxiliary electronic variable momentum
        if "p_aux_var" in saver.keywords:  # and "p_aux_var" in saver.np_data.keys():
            saver.add_dataset("p_aux_var", (_nsteps, 1), "R")

        # Trajectory-resolved auxiliary electronic variable force
        if "f_aux_var" in saver.keywords:  # and "p_aux_var" in saver.np_data.keys():
            saver.add_dataset("f_aux_var", (_nsteps, 1), "R")

    if output_level >= 3:

        # Average adiabatic SH populations (dynamically-consistent)
        if "SH_pop" in saver.keywords:  # and "SH_pop" in saver.np_data.keys():
            saver.add_dataset("SH_pop", (_nsteps, _nadi, 1), "R")

        # Average adiabatic SH populations (raw)
        if "SH_pop_raw" in saver.keywords:  # and "SH_pop_raw" in saver.np_data.keys():
            saver.add_dataset("SH_pop_raw", (_nsteps, _nadi, 1), "R")

        # Average adiabatic density matrices (dynamically-consistent)
        if "D_adi" in saver.keywords:  # and "D_adi" in saver.np_data.keys():
            saver.add_dataset("D_adi", (_nsteps, _nadi, _nadi), "C")

        # Average adiabatic density matrices (raw)
        if "D_adi_raw" in saver.keywords:  # and "D_adi_raw" in saver.np_data.keys():
            saver.add_dataset("D_adi_raw", (_nsteps, _nadi, _nadi), "C")

        # Average diabatic density matrices (dynamically-consistent)
        if "D_dia" in saver.keywords:  # and "D_dia" in saver.np_data.keys():
            saver.add_dataset("D_dia", (_nsteps, _ndia, _ndia), "C")

        # Average diabatic density matrices (raw)
        if "D_dia_raw" in saver.keywords:  # and "D_dia_raw" in saver.np_data.keys():
            saver.add_dataset("D_dia_raw", (_nsteps, _ndia, _ndia), "C")

        # Average coherence indicator matrices
        if "coherence_adi" in saver.keywords:  # and "D_adi" in saver.np_data.keys():
            saver.add_dataset("coherence_adi", (_nsteps, _nadi, _nadi), "R")

        if "coherence_dia" in saver.keywords:  # and "D_adi" in saver.np_data.keys():
            saver.add_dataset("coherence_dia", (_nsteps, _ndia, _ndia), "R")

        # Trajectory-resolved coordinates
        if "q" in saver.keywords:  # and "q" in saver.np_data.keys():
            saver.add_dataset("q", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved momenta
        if "p" in saver.keywords:  # and "p" in saver.np_data.keys():
            saver.add_dataset("p", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved active forces
        if "f" in saver.keywords:  # and "f" in saver.np_data.keys():
            saver.add_dataset("f", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved adiabatic TD-SE amplitudes
        if "Cadi" in saver.keywords:  # and "Cadi" in saver.np_data.keys():
            saver.add_dataset("Cadi", (_nsteps, _ntraj, _nadi), "C")

        # Trajectory-resolved diabatic TD-SE amplitudes
        if "Cdia" in saver.keywords:  # and "Cdia" in saver.np_data.keys():
            saver.add_dataset("Cdia", (_nsteps, _ntraj, _ndia), "C")

        # Trajectory-resolved adiabatic MMST electronic "coordinates"
        if "q_mm" in saver.keywords:  # and "q_mm" in saver.np_data.keys():
            saver.add_dataset("q_mm", (_nsteps, _ntraj, _nadi), "R")

        # Trajectory-resolved adiabatic MMST electronic "momenta"
        if "p_mm" in saver.keywords:  # and "p_mm" in saver.np_data.keys():
            saver.add_dataset("p_mm", (_nsteps, _ntraj, _nadi), "R")

        # Trajectory-resolved quantum momenta
        if "wp_width" in saver.keywords:  # and "p_quant" in saver.np_data.keys():
            saver.add_dataset("wp_width", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved quantum momenta
        if "p_quant" in saver.keywords:  # and "p_quant" in saver.np_data.keys():
            saver.add_dataset("p_quant", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved exact vector potential
        if "VP" in saver.keywords:  # and "p_quant" in saver.np_data.keys():
            saver.add_dataset("VP", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved decoherence forces based on XF
        if "f_xf" in saver.keywords:  # and "f_xf" in saver.np_data.keys():
            saver.add_dataset("f_xf", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved nonclassical forces in QTSH
        if "qtsh_f_nc" in saver.keywords:  # and "f_xf" in saver.np_data.keys():
            saver.add_dataset("qtsh_f_nc", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory-resolved decoherence rates
        if "ave_decoherence_rates" in saver.keywords:  # decoherence time:
            saver.add_dataset("ave_decoherence_rates", (_nsteps, _nadi, _nadi), "R") 

    if output_level >= 4:

        # Trajectory-resolved vibronic Hamiltoninans in the adiabatic representation
        if "hvib_adi" in saver.keywords:  # and "hvib_adi" in saver.np_data.keys():
            saver.add_dataset("hvib_adi", (_nsteps, _ntraj, _nadi, _nadi), "C")

        # Trajectory-resolved vibronic Hamiltoninans in the diabatic representation
        if "hvib_dia" in saver.keywords:  # and "hvib_dia" in saver.np_data.keys():
            saver.add_dataset("hvib_dia", (_nsteps, _ntraj, _ndia, _ndia), "C")

        # Trajectory-resolved time-overlaps of the adiabatic states
        if "St" in saver.keywords:  # and "St" in saver.np_data.keys():
            saver.add_dataset("St", (_nsteps, _ntraj, _nadi, _nadi), "C")

        # Trajectory-resolved diabatic-to-adiabatic transformation matrices
        if "basis_transform" in saver.keywords:  # and "basis_transform" in saver.np_data.keys():
            saver.add_dataset("basis_transform", (_nsteps, _ntraj, _ndia, _nadi), "C")

        # Trajectory-resolved projector matrices (from the raw adiabatic to consistent adiabatic)
        if "projector" in saver.keywords:  # and "projector" in saver.np_data.keys():
            saver.add_dataset("projector", (_nsteps, _ntraj, _nadi, _nadi), "C")

        # Trajectory-resolved auxiliary coordinates
        if "q_aux" in saver.keywords:  # and "hvib_adi" in saver.np_data.keys():
            saver.add_dataset("q_aux", (_nsteps, _ntraj, _nadi, _ndof), "R")

        # Trajectory-resolved auxiliary momenta
        if "p_aux" in saver.keywords:  # and "hvib_adi" in saver.np_data.keys():
            saver.add_dataset("p_aux", (_nsteps, _ntraj, _nadi, _ndof), "R")

        # Trajectory-resolved nabla_phase
        if "nab_phase" in saver.keywords:  # and "hvib_adi" in saver.np_data.keys():
            saver.add_dataset("nab_phase", (_nsteps, _ntraj, _nadi, _ndof), "R")

    if output_level >= 5:
        # Trajectory-resolved derivative coupling vectors
        if "dc1_adi" in saver.keywords:  # and "hvib_adi" in saver.np_data.keys():
            saver.add_dataset("dc1_adi", (_nsteps, _ntraj, _ndof, _nadi, _nadi), "R")


def init_tsh_savers(params, model_params, nsteps, ntraj, nnucl, nadi, ndia):

    # ================ Create savers ==================
    if params["txt_output_level"] > 0 or params["hdf5_output_level"] > 0 or params["mem_output_level"] > 0:

        prefix = params["prefix"]

        # Create an output directory, if not present
        if not os.path.isdir(prefix):
            os.mkdir(prefix)

        # Simulation parameters
        f = open(F"{prefix}/_dyn_params.txt", "w")
        f.write(str(params))
        f.close()

        f = open(F"{prefix}/_model_params.txt", "w")
        f.write(str(model_params))
        f.close()

    if params["txt2_output_level"] > 0:

        prefix2 = params["prefix2"]

        # Create an output directory, if not present
        if not os.path.isdir(prefix2):
            os.mkdir(prefix2)

        # Simulation parameters
        f = open(F"{prefix2}/_dyn_params.txt", "w")
        f.write(str(params))
        f.close()

        f = open(F"{prefix2}/_model_params.txt", "w")
        f.write(str(model_params))
        f.close()

    properties_to_save = params["properties_to_save"]

    _savers = {"hdf5_saver": None, "txt_saver": None, "mem_saver": None, "txt2_saver": None}

    # ====== HDF5 ========
    hdf5_output_level = params["hdf5_output_level"]

    if hdf5_output_level > 0:
        _savers["hdf5_saver"] = data_savers.hdf5_saver(F"{prefix}/data.hdf", properties_to_save)
        _savers["hdf5_saver"].set_compression_level(params["use_compression"], params["compression_level"])
        init_tsh_data(_savers["hdf5_saver"], hdf5_output_level, nsteps, ntraj, nnucl, nadi, ndia)

    # ====== TXT ========
    txt_output_level = params["txt_output_level"]

    if params["txt_output_level"] > 0:
        _savers["txt_saver"] = data_savers.mem_saver(properties_to_save)
        # _savers["txt_saver"].set_compression_level(params["use_compression"], params["compression_level"])
        init_tsh_data(_savers["txt_saver"], txt_output_level, nsteps, ntraj, nnucl, nadi, ndia)

    # ====== TXT2: No intermediate memory allocation ========
    txt2_output_level = params["txt2_output_level"]

    if params["txt2_output_level"] > 0:
        _savers["txt2_saver"] = data_savers.mem_saver(properties_to_save)

        # Here, nsteps is set to 1 since in this type of saver, we only care about the current values,
        # not all the timesteps - that would be too consuming
        init_tsh_data(_savers["txt2_saver"], txt2_output_level, 1, ntraj, nnucl, nadi, ndia)

    # ====== MEM =========
    mem_output_level = params["mem_output_level"]

    if mem_output_level > 0:
        _savers["mem_saver"] = data_savers.mem_saver(properties_to_save)
        init_tsh_data(_savers["mem_saver"], mem_output_level, nsteps, ntraj, nnucl, nadi, ndia)

    return _savers


def save_hdf5_1D(saver, i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type == 0:
        t = i

    # Timestep
    saver.save_scalar(t, "timestep", i)

    # Actual time
    saver.save_scalar(t, "time", dt * i)

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
    if txt_type == 0:
        t = i

    # Timestep
    saver.save_scalar(t, "timestep", i)

    # Actual time
    saver.save_scalar(t, "time", dt * i)

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

    # TC-NBRA average kinetic energy
    if "tcnbra_ekin" in saver.keywords and "tcnbra_ekin" in saver.np_data.keys():
        tcnbra_ekin = dyn_var.compute_tcnbra_ekin()
        saver.save_scalar(t, "tcnbra_ekin", tcnbra_ekin)

    # TC-NBRA average energy of thermostat
    if "tcnbra_thermostat_energy" in saver.keywords and "tcnbra_thermostat_energy" in saver.np_data.keys():
        tcnbra_thermostat_energy = dyn_var.compute_tcnbra_thermostat_energy()
        saver.save_scalar(t, "tcnbra_thermostat_energy", tcnbra_thermostat_energy)

    # KC-RPMD auxiliary electronic variable kinetic energy
    if "ekin_aux_var" in saver.keywords and "ekin_aux_var" in saver.np_data.keys():
        ekin_aux_var = dyn_var.compute_kcrpmd_ekin()
        saver.save_scalar(t, "ekin_aux_var", ekin_aux_var)


def save_hdf5_2D(saver, i, states, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type == 0:
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
    if txt_type == 0:
        t = i

    if "states" in saver.keywords and "states" in saver.np_data.keys():
        # Trajectory-resolved instantaneous adiabatic states
        # Format: saver.add_dataset("states", (_nsteps, _ntraj), "I")
        ntraj = dyn_var.ntraj
        for itraj in range(ntraj):
            saver.save_multi_scalar(t, itraj, "states", dyn_var.act_states[itraj])

    if "states_dia" in saver.keywords and "states_dia" in saver.np_data.keys():
        # Trajectory-resolved instantaneous diabatic states
        # Format: saver.add_dataset("states", (_nsteps, _ntraj), "I")
        ntraj = dyn_var.ntraj
        for itraj in range(ntraj):
            saver.save_multi_scalar(t, itraj, "states_dia", dyn_var.act_states_dia[itraj])

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

    if "sh_pop_dia" in saver.keywords and "sh_pop_dia" in saver.np_data.keys():
        # Average diabatic SH populations
        # Format: saver.add_dataset("sh_pop_dia", (_nsteps, _ndia), "R")
        pops_sh0 = dyn_var.compute_average_sh_pop(0)
        ndia = dyn_var.ndia
        for ist in range(ndia):
            saver.save_multi_scalar(t, ist, "sh_pop_dia", pops_sh0[ist])

    if "sh_pop_adi" in saver.keywords and "sh_pop_adi" in saver.np_data.keys():
        # Average adiabatic SH populations
        # Format: saver.add_dataset("sh_pop_adi", (_nsteps, _nadi), "R")
        pops_sh1 = dyn_var.compute_average_sh_pop(1)
        nadi = dyn_var.nadi
        for ist in range(nadi):
            saver.save_multi_scalar(t, ist, "sh_pop_adi", pops_sh1[ist])

    if "sh_pop_dia_TR" in saver.keywords and "sh_pop_dia_TR" in saver.np_data.keys():
        # Average diabatic SH populations
        # Format: saver.add_dataset("sh_pop_dia", (_nsteps, _ndia), "R")
        pops_sh0 = dyn_var.compute_average_sh_pop_TR(0)
        ndia = dyn_var.ndia
        for ist in range(ndia):
            saver.save_multi_scalar(t, ist, "sh_pop_dia_TR", pops_sh0[ist])

    if "sh_pop_adi_TR" in saver.keywords and "sh_pop_adi_TR" in saver.np_data.keys():
        # Average adiabatic SH populations
        # Format: saver.add_dataset("sh_pop_adi", (_nsteps, _nadi), "R")
        pops_sh1 = dyn_var.compute_average_sh_pop_TR(1)
        nadi = dyn_var.nadi
        for ist in range(nadi):
            saver.save_multi_scalar(t, ist, "sh_pop_adi_TR", pops_sh1[ist])

    if "mash_pop_dia" in saver.keywords and "mash_pop_dia" in saver.np_data.keys():
        # Average diabatic MASH populations
        # Format: saver.add_dataset("mash_pop_dia", (_nsteps, _ndia), "R")
        pops_sh0 = dyn_var.compute_average_mash_pop(0)
        ndia = dyn_var.ndia
        for ist in range(ndia):
            saver.save_multi_scalar(t, ist, "mash_pop_dia", pops_sh0[ist])

    if "mash_pop_adi" in saver.keywords and "mash_pop_adi" in saver.np_data.keys():
        # Average adiabatic MASH populations
        # Format: saver.add_dataset("mash_pop_adi", (_nsteps, _nadi), "R")
        pops_sh1 = dyn_var.compute_average_mash_pop(1)
        nadi = dyn_var.nadi
        for ist in range(nadi):
            saver.save_multi_scalar(t, ist, "mash_pop_adi", pops_sh1[ist])

    if "fssh3_average_errors" in saver.keywords and "fssh3_average_errors" in saver.np_data.keys():
        # Average errors in FSSH3
        # Format: saver.add_dataset("fssh3_average_errors", (_nsteps, 5), "R")
        fssh3_average_errors = dyn_var.get_fssh3_average_errors()
        for k in range(5):
            saver.save_multi_scalar(t, k, "fssh3_average_errors", fssh3_average_errors[k])

    if "y_aux_var" in saver.keywords and "y_aux_var" in saver.np_data.keys():
        # Trajectory-resolved auxiliary electronic variable coordinate
        # Format: saver.add_dataset("y_aux_var", (_nsteps, 1), "R")
        y_aux_var = dyn_var.get_y_aux_var()
        saver.save_multi_scalar(t, 0, "y_aux_var", y_aux_var[0])

    if "p_aux_var" in saver.keywords and "p_aux_var" in saver.np_data.keys():
        # Trajectory-resolved auxiliary electronic variable coordinate
        # Format: saver.add_dataset("p_aux_var", (_nsteps, 1), "R")
        p_aux_var = dyn_var.get_p_aux_var()
        saver.save_multi_scalar(t, 0, "p_aux_var", p_aux_var[0])

    if "f_aux_var" in saver.keywords and "f_aux_var" in saver.np_data.keys():
        # Trajectory-resolved auxiliary electronic variable coordinate
        # Format: saver.add_dataset("f_aux_var", (_nsteps, 1), "R")
        f_aux_var = dyn_var.get_f_aux_var()
        saver.save_multi_scalar(t, 0, "f_aux_var", f_aux_var[0])


def save_hdf5_3D(saver, i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type == 0:
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
    if txt_type == 0:
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

    # Trajectory-resolved MMST electronic "coordinate"
    # Format: saver.add_dataset("q_mm", (_nsteps, _ntraj, _nadi), "C")
    if "q_mm" in saver.keywords and "q_mm" in saver.np_data.keys():
        q_mm = dyn_var.get_q_mm()
        saver.save_matrix(t, "q_mm", q_mm.T())

    # Trajectory-resolved MMST electronic "momentum"
    # Format: saver.add_dataset("p_mm", (_nsteps, _ntraj, _nadi), "C")
    if "p_mm" in saver.keywords and "p_mm" in saver.np_data.keys():
        p_mm = dyn_var.get_p_mm()
        saver.save_matrix(t, "p_mm", p_mm.T())

    # Average adiabatic coherence indicator matrices
    # Format: saver.add_dataset("coherence_adi", (_nsteps, _nadi, _nadi), "R")
    if "coherence_adi" in saver.keywords and "coherence_adi" in saver.np_data.keys():
        coh_adi = dyn_var.compute_coherence_indicator(1)
        saver.save_matrix(t, "coherence_adi", coh_adi)

    # Average diabatic coherence indicator matrices
    # Format: saver.add_dataset("coherence_dia", (_nsteps, _ndia, _ndia), "R")
    if "coherence_dia" in saver.keywords and "coherence_dia" in saver.np_data.keys():
        coh_dia = dyn_var.compute_coherence_indicator(0)
        saver.save_matrix(t, "coherence_dia", coh_dia)

    # Wavepacket width
    # Format: saver.add_dataset("wp_width", (_nsteps, _ntraj, _dof), "R")
    if "wp_width" in saver.keywords and "wp_width" in saver.np_data.keys():
        wp_width = dyn_var.get_wp_width()
        saver.save_matrix(t, "wp_width", wp_width.T())

    # Trajectory-resolved quantum momenta
    # Format: saver.add_dataset("p_quant", (_nsteps, _ntraj, _dof), "R")
    if "p_quant" in saver.keywords and "p_quant" in saver.np_data.keys():
        p_quant = dyn_var.get_p_quant()
        saver.save_matrix(t, "p_quant", p_quant.T())

    # Trajectory-resolved exact vector potential
    # Format: saver.add_dataset("VP", (_nsteps, _ntraj, _dof), "R")
    if "VP" in saver.keywords and "VP" in saver.np_data.keys():
        VP = dyn_var.get_VP()
        saver.save_matrix(t, "VP", VP.T())

    # Trajectory-resolved XF-based decoherence force
    # Format: saver.add_dataset("f_xf", (_nsteps, _ntraj, _dof), "R")
    if "f_xf" in saver.keywords and "f_xf" in saver.np_data.keys():
        f_xf = dyn_var.get_f_xf()
        saver.save_matrix(t, "f_xf", f_xf.T())

    # Nonclassical force in QTSH
    # Format: saver.add_dataset("qtsh_f_nc", (_nsteps, _ntraj, _dof), "R")
    if "qtsh_f_nc" in saver.keywords and "qtsh_f_nc" in saver.np_data.keys():
        qtsh_f_nc = dyn_var.get_qtsh_f_nc()
        saver.save_matrix(t, "qtsh_f_nc", qtsh_f_nc.T())

    # Average decoherence rate
    # Format: saver.add_dataset("ave_decoherence_rates", (_nsteps, _nadi, _nadi), "R") 
    if "ave_decoherence_rates" in saver.keywords and "ave_decoherence_rates" in saver.np_data.keys():
        ave_decoherence_rates = dyn_var.get_ave_decoherence_rates()
        saver.save_matrix(t, "ave_decoherence_rates", ave_decoherence_rates.T())


def save_hdf5_4D(
        saver,
        i,
        tr,
        hvib_adi,
        hvib_dia,
        St,
        U,
        projector,
        q_aux=None,
        p_aux=None,
        nab_phase=None,
        txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type == 0:
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

    # Trajectory-resolved auxiliary coordinates
    # Format: saver.add_dataset("q_aux", (_nsteps, _ntraj, _nadi, _ndof), "R")
    if "q_aux" in saver.keywords and "q_aux" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "q_aux", q_aux)

    # Trajectory-resolved auxiliary momenta
    # Format: saver.add_dataset("p_aux", (_nsteps, _ntraj, _nadi, _ndof), "R")
    if "p_aux" in saver.keywords and "p_aux" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "p_aux", p_aux)

    # Trajectory-resolved spatial derivatives of coefficient phase
    # Format: saver.add_dataset("nab_phase", (_nsteps, _ntraj, _nadi, _ndof), "R")
    if "nab_phase" in saver.keywords and "nab_phase" in saver.np_data.keys():
        saver.save_multi_matrix(t, tr, "nab_phase", nab_phase)


def save_hdf5_5D(saver, i, tr, idof, dc1_adi, txt_type=0):
    """
    saver - can be either hdf5_saver or mem_saver

    txt_type ( int ): 0 - standard, all the timesteps, 1 - only the current one

    """

    t = 0
    if txt_type == 0:
        t = i

    # Trajectory-resolved vibronic Hamiltoninans in the adiabatic representation
    # Format: saver.add_dataset("dc1_adi", (_nsteps, _ntraj, _nadi, _nadi, _ndof), "R")
    if "dc1_adi" in saver.keywords and "dc1_adi" in saver.np_data.keys():
        saver.save_multi_matrix2(t, tr, idof, "dc1_adi", dc1_adi)


def save_tsh_data_123(_savers, params,
                      i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, states,
                      pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia
                      ):

    hdf5_output_level = params["hdf5_output_level"]
    mem_output_level = params["mem_output_level"]
    txt_output_level = params["txt_output_level"]
    txt2_output_level = params["txt2_output_level"]

    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"] * nsteps)

    if i % print_freq == 0:
        print(F" step= {i}")

    if hdf5_output_level >= 1 and _savers["hdf5_saver"] is not None:
        save_hdf5_1D(_savers["hdf5_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)

    if mem_output_level >= 1 and _savers["mem_saver"] is not None:
        save_hdf5_1D(_savers["mem_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)

    if txt_output_level >= 1 and _savers["txt_saver"] is not None:
        save_hdf5_1D(_savers["txt_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)

    if txt2_output_level >= 1 and _savers["txt2_saver"] is not None:
        save_hdf5_1D(_savers["txt2_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, 1)

    if hdf5_output_level >= 2 and _savers["hdf5_saver"] is not None:
        save_hdf5_2D(_savers["hdf5_saver"], i, states)

    if mem_output_level >= 2 and _savers["mem_saver"] is not None:
        save_hdf5_2D(_savers["mem_saver"], i, states)

    if txt_output_level >= 2 and _savers["txt_saver"] is not None:
        save_hdf5_2D(_savers["txt_saver"], i, states)

    if txt2_output_level >= 2 and _savers["txt2_saver"] is not None:
        save_hdf5_2D(_savers["txt2_saver"], i, states, 1)

    if hdf5_output_level >= 3 and _savers["hdf5_saver"] is not None:
        save_hdf5_3D(_savers["hdf5_saver"], i, pops, pops_raw, dm_adi,
                     dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)

    if mem_output_level >= 3 and _savers["mem_saver"] is not None:
        save_hdf5_3D(_savers["mem_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)

    if txt_output_level >= 3 and _savers["txt_saver"] is not None:
        save_hdf5_3D(_savers["txt_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)

    if txt2_output_level >= 3 and _savers["txt2_saver"] is not None:
        save_hdf5_3D(_savers["txt2_saver"], i, pops, pops_raw, dm_adi,
                     dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia, 1)


def save_tsh_data_1234_new(_savers, params, i, dyn_var, ham):

    hdf5_output_level = params["hdf5_output_level"]
    mem_output_level = params["mem_output_level"]
    txt_output_level = params["txt_output_level"]
    txt2_output_level = params["txt2_output_level"]

    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"] * nsteps)

    if i % print_freq == 0:
        print(F" step= {i}")

    # =========== Using: def save_hdf5_1D_new(saver, i, params, dyn_var, ham, txt_type=0) =========

    if hdf5_output_level >= 1 and _savers["hdf5_saver"] is not None:
        save_hdf5_1D_new(_savers["hdf5_saver"], i, params, dyn_var, ham)

    if mem_output_level >= 1 and _savers["mem_saver"] is not None:
        # save_hdf5_1D(_savers["mem_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)
        save_hdf5_1D_new(_savers["mem_saver"], i, params, dyn_var, ham)

    if txt_output_level >= 1 and _savers["txt_saver"] is not None:
        # save_hdf5_1D(_savers["txt_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC)
        save_hdf5_1D_new(_savers["txt_saver"], i, params, dyn_var, ham)

    if txt2_output_level >= 1 and _savers["txt2_saver"] is not None:
        # save_hdf5_1D(_savers["txt2_saver"], i, dt, Ekin, Epot, Etot, dEkin, dEpot, dEtot, Etherm, E_NHC, 1)
        save_hdf5_1D_new(_savers["txt2_saver"], i, params, dyn_var, ham, 1)

    # =========== Using : def save_hdf5_2D_new(saver, i, dyn_var, txt_type=0) =====================
    if hdf5_output_level >= 2 and _savers["hdf5_saver"] is not None:
        # save_hdf5_2D(_savers["hdf5_saver"], i, states)
        save_hdf5_2D_new(_savers["hdf5_saver"], i, dyn_var, ham)

    if mem_output_level >= 2 and _savers["mem_saver"] is not None:
        # save_hdf5_2D(_savers["mem_saver"], i, states)
        save_hdf5_2D_new(_savers["mem_saver"], i, dyn_var, ham)

    if txt_output_level >= 2 and _savers["txt_saver"] is not None:
        # save_hdf5_2D(_savers["txt_saver"], i, states)
        save_hdf5_2D_new(_savers["txt_saver"], i, dyn_var, ham)

    if txt2_output_level >= 2 and _savers["txt2_saver"] is not None:
        # save_hdf5_2D(_savers["txt2_saver"], i, states, 1)
        save_hdf5_2D_new(_savers["txt2_saver"], i, dyn_var, ham, 1)

    # ============ Using: def save_hdf5_3D_new(saver, i, dyn_var, txt_type=0) ==========================

    if hdf5_output_level >= 3 and _savers["hdf5_saver"] is not None:
        # save_hdf5_3D(_savers["hdf5_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)
        save_hdf5_3D_new(_savers["hdf5_saver"], i, dyn_var)

    if mem_output_level >= 3 and _savers["mem_saver"] is not None:
        # save_hdf5_3D(_savers["mem_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)
        save_hdf5_3D_new(_savers["mem_saver"], i, dyn_var)

    if txt_output_level >= 3 and _savers["txt_saver"] is not None:
        # save_hdf5_3D(_savers["txt_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia)
        save_hdf5_3D_new(_savers["txt_saver"], i, dyn_var)

    if txt2_output_level >= 3 and _savers["txt2_saver"] is not None:
        # save_hdf5_3D(_savers["txt2_saver"], i, pops, pops_raw, dm_adi, dm_adi_raw, dm_dia, dm_dia_raw, q, p, Cadi, Cdia, 1)
        save_hdf5_3D_new(_savers["txt2_saver"], i, dyn_var, 1)

    # ============ Using: def save_hdf5_4D_new(saver, i, tr, dyn_var, ham, txt_type=0) ====================

    is_nbra = params["is_nbra"]
    time_overlap_method = params["time_overlap_method"]

    tr_range = []
    if is_nbra == 1:
        tr_range = [0]
    else:
        tr_range = range(dyn_var.ntraj)

    for tr in tr_range:

        U = ham.get_basis_transform(Py2Cpp_int([0, tr]))
        St = ham.get_time_overlap_adi(Py2Cpp_int([0, tr]))
        T = dyn_var.get_proj_adi(tr)

        q_aux, p_aux, nab_phase = None, None, None
        if dyn_var.shxf_vars_status == 1 or dyn_var.mqcxf_vars_status == 1:
            q_aux = dyn_var.get_coords_aux(tr)
            p_aux = dyn_var.get_momenta_aux(tr)
            nab_phase = dyn_var.get_nab_phase(tr)

        if hdf5_output_level >= 4 and _savers["hdf5_saver"] is not None:
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr]))
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr]))
            save_hdf5_4D(_savers["hdf5_saver"], i, tr, hvib_adi, hvib_dia, St, U, T, q_aux, p_aux, nab_phase)

        if mem_output_level >= 4 and _savers["mem_saver"] is not None:
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr]))
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr]))
            save_hdf5_4D(_savers["mem_saver"], i, tr, hvib_adi, hvib_dia, St, U, T, q_aux, p_aux, nab_phase)

        if txt_output_level >= 4 and _savers["txt_saver"] is not None:
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr]))
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr]))
            save_hdf5_4D(_savers["txt_saver"], i, tr, hvib_adi, hvib_dia, St, U, T, q_aux, p_aux, nab_phase)

        if txt2_output_level >= 4 and _savers["txt2_saver"] is not None:
            hvib_adi = ham.get_hvib_adi(Py2Cpp_int([0, tr]))
            hvib_dia = ham.get_hvib_dia(Py2Cpp_int([0, tr]))
            save_hdf5_4D(_savers["txt2_saver"], i, tr, hvib_adi, hvib_dia, St, U, T, q_aux, p_aux, nab_phase, 1)

    # ============= Using save_hdf5_5D(saver, i, tr, idof, dc1) =========================
    for tr in tr_range:
        for idof in range(dyn_var.ndof):

            dc1 = ham.get_dc1_adi(idof, (Py2Cpp_int([0, tr])))  # .real()
            dc1 = dc1.real()
            # T = dyn_var.get_proj_adi(tr)
            # dc1 = (T.H() * dc1 * T).real()

            if hdf5_output_level >= 5 and _savers["hdf5_saver"] is not None:
                save_hdf5_5D(_savers["hdf5_saver"], i, tr, idof, dc1)

            if mem_output_level >= 5 and _savers["mem_saver"] is not None:
                save_hdf5_5D(_savers["mem_saver"], i, tr, idof, dc1)

            if txt_output_level >= 5 and _savers["txt_saver"] is not None:
                save_hdf5_5D(_savers["txt_saver"], i, tr, idof, dc1)

            if txt2_output_level >= 5 and _savers["txt2_saver"] is not None:
                save_hdf5_5D(_savers["txt_saver"], i, tr, idof, dc1)


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
        add_doublelist2file("%s/energies.txt" % (prefix), i * dt, [Ekin, Epot, Etot, dEkin, dEpot, dEtot])
        add_cmatrix2file("%s/D_adi.txt" % (prefix), i * dt, dm_adi)
        add_cmatrix2file("%s/D_dia.txt" % (prefix), i * dt, dm_dia)
        add_matrix2file("%s/SH_pop.txt" % (prefix), i * dt, pops)

    # File output
    if file_output_level >= 2:
        add_matrix2file("%s/q.txt" % (prefix), i * dt, q)
        add_matrix2file("%s/p.txt" % (prefix), i * dt, p)
        add_cmatrix2file("%s/C_adi.txt" % (prefix), i * dt, Cadi)
        add_cmatrix2file("%s/C_dia.txt" % (prefix), i * dt, Cdia)
        add_intlist2file("%s/states.txt" % (prefix), i * dt, states)


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
        add_cmatrix2file("%s/Hvib_adi_%i.txt" % (prefix, tr), i * dt, hvib_adi)
        add_cmatrix2file("%s/Hvib_dia_%i.txt" % (prefix, tr), i * dt, hvib_dia)
        add_cmatrix2file("%s/St_%i.txt" % (prefix, tr), i * dt, St)
        add_cmatrix2file("%s/basis_transform_%i.txt" % (prefix, tr), i * dt, U)


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
        res = data_read.get_data_from_file2("%s/energies.txt" % (prefix), range(0, 7))
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
