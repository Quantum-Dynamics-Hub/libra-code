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

    if output_level>=3:

        # Trajectory positions
        if "traj_q" in saver.keywords: # and "traj_q" in saver.np_data.keys():
            saver.add_dataset("traj_q", (_nsteps, _ndof, _ntraj), "R")

        # Trajectory x-dependent phases
        if "traj_p" in saver.keywords: # and "traj_p" in saver.np_data.keys():
            saver.add_dataset("traj_p", (_nsteps, _ndof, _ntraj), "R")

        # Trajectory widths
        if "traj_a" in saver.keywords: # and "traj_a" in saver.np_data.keys():
            saver.add_dataset("traj_a", (_nsteps, _ndof, _ntraj), "R")

        # Trajectory x-independent phases
        if "traj_s" in saver.keywords: # and "traj_s" in saver.np_data.keys():
            saver.add_dataset("traj_s", (_nsteps, _ndof, _ntraj), "R")


    if output_level>=4:

        # Trajectory basis coefficients
        if "coeffs" in saver.keywords: # and "coeffs" in saver.np_data.keys():
            saver.add_dataset("coeffs", (_nsteps, _ntraj), "C")

def init_qtag_savers(params, model_params, nsteps, ntraj, ndof, nstates):
    #====== CREATE DIRECTORIES ========
    if params["hdf5_output_level"] > 0 or params["mem_output_level"] > 0:

        prefix = params["prefix"]

        # Create an output directory, if not present
        if not os.path.isdir(prefix):
            os.mkdir(prefix)

        # Simulation parameters
        f = open(F"{prefix}/_dyn_params.txt","w")
        f.write( str(params) );  f.close()

        f = open(F"{prefix}/_model_params.txt","w")
        f.write( str(model_params) );  f.close()

    #====== HDF5 ========
    _savers = {"hdf5_saver":None, "mem_saver":None}

    hdf5_output_level = params["hdf5_output_level"]
    properties_to_save = params["properties_to_save"]

    if hdf5_output_level > 0:
        _savers["hdf5_saver"] = data_savers.hdf5_saver(F"{prefix}/data.hdf", properties_to_save)
        _savers["hdf5_saver"].set_compression_level(params["use_compression"], params["compression_level"])
        init_qtag_data(_savers["hdf5_saver"], hdf5_output_level, nsteps, ntraj, ndof, nstates)

    return(_savers)

def save_qtag_hdf5_lvl1(saver,params,step,Etot,dEtot):

    dt = params['dt']

    #Time
    saver.save_scalar(step, "time", dt*step)

    #Energy
    saver.save_scalar(step, "Etot", Etot)

    #Energy difference from t=0
    saver.save_scalar(step, "dEtot", dEtot)

def save_qtag_hdf5_lvl2(saver,params,step,pops):

    #Surface populations
    saver.save_matrix(step, "pops", pops)

def save_qtag_hdf5_lvl3(saver,params,step,qvals,pvals,avals,svals):

    #Trajectory positions
    saver.save_matrix(step, "traj_q", qvals)

    #Trajectory momenta
    saver.save_matrix(step, "traj_p", pvals)

    #Trajectory widths
    saver.save_matrix(step, "traj_a", avals)

    #Trajectory phases
    saver.save_matrix(step, "traj_s", svals)

def save_qtag_hdf5_lvl4(saver,params,step,ctot):

    #Trajectory coefficients
    saver.save_matrix(step, "coeffs", ctot.T())

