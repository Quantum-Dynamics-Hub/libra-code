#*********************************************************************************
#* Copyright (C) 2021-2022 Matthew Dutra, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#***********************************************************************************
"""
..module:: save
  :platform: Unix, Windows
  :synopsis: This module contains functions for saving variables of QTAG workflow

..moduleauthors :: Matthew Dutra, Alexey Akimov
"""

import os
import sys
import math
import copy

from liblibra_core import *

import util.libutil as comn
import libra_py.data_savers as data_savers

def init_qtag_data(saver, output_level, _nsteps, _ntraj, _ndof, _nstates):
    """
    saver - can be either hdf5_saver or mem_saver

    """

    if output_level>=1:

        # Time axis
        saver.add_dataset("time", (_nsteps,) , "R")

        # Kinetic energy
#        saver.add_dataset("Ekin", (_nsteps,) , "R")

        # Potential energy
#        saver.add_dataset("Epot", (_nsteps,) , "R")

        # Total energy
        saver.add_dataset("Etot", (_nsteps,) , "R")

        # Fluctuation of kinetic energy
#        saver.add_dataset("dEkin", (_nsteps,) , "R")

        # Fluctuation of potential energy
#        saver.add_dataset("dEpot", (_nsteps,) , "R")

        # Fluctuation of total energy
        saver.add_dataset("dEtot", (_nsteps,) , "R")

    if output_level>=2:

        # State populations
        if "pops" in saver.keywords: # and "pops" in saver.np_data.keys():
            saver.add_dataset("pops", (_nsteps, _nstates), "R")

        # Trajectory basis coefficients
        if "coeffs" in saver.keywords: # and "coeffs" in saver.np_data.keys():
            saver.add_dataset("coeffs", (_nsteps, 1, _ntraj), "C")

    if output_level>=3:

        # Trajectory positions
        if "q" in saver.keywords: # and "q" in saver.np_data.keys():
            saver.add_dataset("q", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory x-dependent phases
        if "p" in saver.keywords: # and "p" in saver.np_data.keys():
            saver.add_dataset("p", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory widths
        if "a" in saver.keywords: # and "a" in saver.np_data.keys():
            saver.add_dataset("a", (_nsteps, _ntraj, _ndof), "R")

        # Trajectory x-independent phases
        if "s" in saver.keywords: # and "s" in saver.np_data.keys():
            saver.add_dataset("s", (_nsteps, _ntraj, _ndof), "R")

def init_qtag_savers(params, model_params, nsteps, ntraj, ndof, nstates):
    #====== CREATE DIRECTORIES ========
    if params["txt2_output_level"] > 0 or params["hdf5_output_level"] > 0 or params["mem_output_level"] > 0:
        
        prefix = params["prefix"]
        properties_to_save = params["properties_to_save"]

        # Create an output directory, if not present
        if not os.path.isdir(prefix):
            os.mkdir(prefix)
        else:
            if os.path.isdir(prefix+"/wfc"):
                subdir=prefix+"/wfc"
                for file in os.listdir(subdir):
                    os.remove(os.path.join(subdir,file))
                os.rmdir(subdir)

            for file in os.listdir(prefix):
                os.remove(os.path.join(prefix,file))

        # Simulation parameters
        f = open(F"{prefix}/_dyn_params.txt","w")
        f.write( str(params) );  f.close()

        f = open(F"{prefix}/_model_params.txt","w")
        f.write( str(model_params) );  f.close()

        _savers = {"txt2_saver":None, "hdf5_saver":None, "mem_saver":None}

    #====== txt2 ========

    txt2_output_level = params["txt2_output_level"]

    if params["txt2_output_level"] > 0:
        _savers["txt2_saver"] = data_savers.mem_saver(properties_to_save)

        # Here, nsteps is set to 1 since in this type of saver, we only care about the current values,
        # not all the timesteps - that would be too consuming
        init_qtag_data(_savers["txt2_saver"], txt2_output_level, 1, ntraj, ndof, nstates)

    #====== HDF5 ========

    hdf5_output_level = params["hdf5_output_level"]

    if hdf5_output_level > 0:
        _savers["hdf5_saver"] = data_savers.hdf5_saver(F"{prefix}/data.hdf", properties_to_save)
        _savers["hdf5_saver"].set_compression_level(params["use_compression"], params["compression_level"])
        init_qtag_data(_savers["hdf5_saver"], hdf5_output_level, nsteps, ntraj, ndof, nstates)

    return(_savers)

def save_qtag_hdf5_1D(saver,dt,step,Etot,dEtot,txt_type=0):


    t = 0
    if txt_type==0:
        t = step

    #Time
    saver.save_scalar(t, "time", dt*step)

    #Energy
    saver.save_scalar(t, "Etot", Etot)

    #Energy difference from t=0
    saver.save_scalar(t, "dEtot", dEtot)

def save_qtag_hdf5_2D(saver,nstates,step,pops,coeffs,txt_type=0):

    t = 0
    if txt_type==0:
        t = step

    #Surface populations
    for state in range(nstates):
        saver.save_multi_scalar(t, state, "pops", pops[state])

    #Trajectory coefficients
    saver.save_matrix(t, "coeffs", coeffs.T())

def save_qtag_hdf5_3D(saver,step,qvals,pvals,avals,svals,txt_type=0):

    t = 0
    if txt_type==0:
        t = step

    #Trajectory positions
    saver.save_matrix(t, "q", qvals.T())

    #Trajectory momenta
    saver.save_matrix(t, "p", pvals.T())

    #Trajectory widths
    saver.save_matrix(t, "a", avals.T())

    #Trajectory phases
    saver.save_matrix(t, "s", svals.T())

def save_qtag_data(_savers, params,
                      step, Etot, dEtot, pops, coeffs, qvals, pvals, avals, svals):

    hdf5_output_level = params["hdf5_output_level"]
    txt2_output_level = params["txt2_output_level"]

    nsteps = params["nsteps"]
    nstates = len(params["states"])
    print_freq = int(params["progress_frequency"]*nsteps)

    #======LEVEL 1======
    if hdf5_output_level>=1 and _savers["hdf5_saver"]!=None:
        save_qtag_hdf5_1D(_savers["hdf5_saver"], params['dt'], step, Etot, dEtot)

    if txt2_output_level>=1 and _savers["txt2_saver"]!=None:
        save_qtag_hdf5_1D(_savers["txt2_saver"], params['dt'], step, Etot, dEtot, 1)

    #======LEVEL 2======
    if hdf5_output_level>=2 and _savers["hdf5_saver"]!=None:
        save_qtag_hdf5_2D(_savers["hdf5_saver"], nstates, step, pops, coeffs)

    if txt2_output_level>=2 and _savers["txt2_saver"]!=None:
        save_qtag_hdf5_2D(_savers["txt2_saver"], nstates, step, pops, coeffs, 1)

    #======LEVEL 3======
    if hdf5_output_level>=3 and _savers["hdf5_saver"]!=None:
        save_qtag_hdf5_3D(_savers["hdf5_saver"], step, qvals, pvals, avals, svals)

    if txt2_output_level>=3 and _savers["txt2_saver"]!=None:
        save_qtag_hdf5_3D(_savers["txt2_saver"], step, qvals, pvals, avals, svals, 1)


