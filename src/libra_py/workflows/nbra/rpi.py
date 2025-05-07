# *********************************************************************************
# * Copyright (C) 2025 Daeho Han and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/

"""
.. module:: rpi
   :platform: Unix, Windows
   :synopsis: this module implements restricted path integral (RPI) calculations

.. moduleauthor:: Daeho Han, Alexey V. Akimov

"""

import os
import sys
import numpy as np
import scipy.sparse as sp
import h5py
import time
import multiprocessing as mp
if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import libra_py.units as units
from libra_py import data_conv
import libra_py.data_read as data_read
import util.libutil as comn
import libra_py.dynamics.tsh.compute as tsh_dynamics
    
def run_patch_dyn_serial(_rpi_params, _dyn_params, compute_model, _model_params, _init_elec, _init_nucl, ibatch, ipatch, istate):
    """
    This function runs a single patch dynamics calculation for the time segment from istep to fstep and a given initial state.
    For conducting a coherent Ehrenfest propagation, `libra_py.dynamics.tsh.compute.generic_recipe` is used.

    Args:
        
        _rpi_params ( dictionary ): parameters controlling the patch dynamics in the RPI calculation
            Can contain:

            * **rpi_params["run_slurm"]** ( bool ): Whether to use the slurm environment to submit the jobs using the submit_template file. 
                If it is set to False, it will run the calculations on the active session but multiple jobs will be run on the current active session.

            * **rpi_params["submit_template"]** ( string ): The path of a template slurm submit file.

            * **rpi_params["submission_exe"]** ( string ): The submission executable

            * **rpi_params["run_python_file"]** ( string ): The path of a template running script.
            
            * **rpi_params["python_exe"]** ( string ): The Python executable

            * **rpi_params["iconds"]** ( list of ints ): The list of initial geometry indices of the trajectory segment used for setting model_params. 
                Each initial condition in this list characterizes each batch.

            * **rpi_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_params["npatches"]** ( int ): The number of patches. The time duration of each patch dynamics will be int(nsteps/npatches) + 1. 
                The additional single step is because the tsh recipe (see libra_py.dynamics.tsh.compute for details) used for the patch dynamics 
                saves the dynamics information before the electron-nuclear propagation in each MD loop.

            * **rpi_params["nstates"]** ( int ): The number of electronic states.

            * **rpi_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_params["prefix"]** ( string ): The prefix of directories having the output patch dynamics. 
                The patch dynamics data are characterized by three numbers - the batch index, the patch index, the initial state index.
                Thus, the full name of output directory is expressed by these three numbers and their limits as follows.
                  rpi_params["prefix"] + F"n{len(rpi_params["iconds"])}" + "_ibatch{X}" + _n{rpi_params["npatches"]} + F"_ipatch{Y}" \
                        + F"_n{rpi_params["nstates"]}" + F"_istate{Z}"
                  X = 0 to len(rpi_params["iconds"]) - 1; Y = 0 to rpi_params["npatches"] - 1; Z = 0 to rpi_params["nstates"] - 1
                An NBRA trajectory from an initial geometry defined in rpi_params['iconds'] is called a batch. 
        
        The rest of parameters are for setting patch dynamics conducted by `libra_py.dynamics.tsh.compute.generic_recipe`.
        
        dyn_params ( dictionary ): control parameters for running the dynamics
            see the documentation of the :func:`libra_py.dynamics.tsh.run_dynamics` function

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations

        _model_params ( dictionary ): contains the selection of a model and the parameters
            for that model Hamiltonian. In addition to all the parameters, should contain the key

            * **_model_params["model0"]** ( int ): the selection of the model to be used to compute diabatic-to-adiabatic
                transformation. The particular value chosen here depends on the way the `compute_model` functions is defined

               :Note: the function selected should be able to generate the transformation matrix.

        _init_elec ( dictionary ): control parameters for initialization of electronic variables
            see the documentation of the :func:`libra_py.dynamics.init_electronic_dyn_var` function

        _init_nucl (dict): controls how to sample nuclear DOFs from what the user provides. Can contain:
          * **_init_nucl["init_type"]** (int) : the type of the sampling:
            - 0 : Coords and momenta are set exactly to the given value
            - 1 : Coords are set, momenta are sampled
            - 2 : Coords are sampled, momenta are set
            - 3 : Both coords and momenta are samples [ default ]
          * **_init_nucl["ndof"]** (int) : the number of nuclear DOFs [default: 1]
          * **_init_nucl["q"]** ( list of `ndof` doubles ): average coordinate of all trajectories, for each nuclear DOF,
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["p"]** ( list of `ndof` doubles ): average momentum of all trajectories, for each nuclear DOF,
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["mass"]** (list of `ndof` doubles): the mass for each nuclear DOF, in a.u. of mass [ default: [2000.0] ]
          * **_init_nucl["force_constant"]** (list of `ndof` doubles) : the force constant of the harmonic oscillator in each
                  nuclear DOF, determins the width of the coordinate and momentum samplings [ default: [0.001] ]

        ibatch ( int ): the initial geometry index of rpi_params["iconds"]
        
        ipatch ( int ): the patch index, running from 0 to rpi_params["npatches"] - 1

        istate ( int ): the initial electronic state for this patch dynamics

    Return:
        None: but performs the action
    """
    rpi_params = dict(_rpi_params)
    model_params = dict(_model_params)
    dyn_params = dict(_dyn_params)
    init_elec = dict(_init_elec)
    init_nucl = dict(_init_nucl)
    
    critical_params = []
    default_params = {"run_slurm": False, "submit_template": 'submit_template.slm', "submission_exe": 'sbatch',
                      "run_python_file": 'run_template.py', "python_exe": 'python',  
                      "nsteps": 1, "npatches": 1, "nstates":2, "iconds": [0], "dt": 1.0*units.fs2au, "prefix": 'out_'}
    
    comn.check_input(rpi_params, default_params, critical_params)
    
    # ========== Model params ==========
    comn.check_input(model_params, {}, ["model0"])

    # ========== Nuclear params ==========
    # First check the inputs for the initialization of nuclear variables
    comn.check_input(init_nucl, {"init_type": 3, "ndof": 1, "q": [0.0], "p": [0.0], "mass": [2000.0],
                                 "force_constant": [0.01], "q_width": [1.0], "p_width": [1.0]}, [])
    ndof = len(init_nucl["mass"])

    # ========== Dynamics params ==========
    # Now, we can use the information of ndof to define the default quantum_dofs
    comn.check_input(dyn_params, {"rep_tdse": 1, "rep_sh": 1, "is_nbra": 0, "direct_init": 0, "ntraj": 1,
                                  "quantum_dofs": list(range(ndof))
                                  }, [])
    
    nsteps, dt = rpi_params["nsteps"], rpi_params["dt"]
    iconds, npatches, nstates = rpi_params["iconds"], rpi_params["npatches"], rpi_params["nstates"]

    # Adjust icond according to the batch and patch indices
    icond = iconds[ibatch] + ipatch * int( nsteps / npatches )
   
    dyn_params.update({"nsteps": int( nsteps / npatches ) + 1, "dt": dt, "icond":icond})

    istates = [0.0]*nstates; istates[istate] = 1.0
    init_elec.update({"istate":istate, "istates": istates})

    rnd = Random()

    res = tsh_dynamics.generic_recipe(dyn_params, compute_model, model_params, init_elec, init_nucl, rnd)

def generic_patch_dyn(_rpi_params, _dyn_params, compute_model, _model_params, _init_elec, _init_nucl):
    """
    This function performs patch dynamics calculations based on the restricted path integral (RPI) approach.

    Args:

        _rpi_params ( dictionary ): parameters controlling the patch dynamics in the RPI calculation
            Can contain:

            * **rpi_params["run_slurm"]** ( bool ): Whether to use the slurm environment to submit the jobs using the submit_template file. 
                If it is set to False, it will run the calculations on the active session but multiple jobs will be run on the current active session.

            * **rpi_params["submit_template"]** ( string ): The path of a template slurm submit file.

            * **rpi_params["submission_exe"]** ( string ): The submission executable

            * **rpi_params["run_python_file"]** ( string ): The path of a template running script.
            
            * **rpi_params["python_exe"]** ( string ): The Python executable

            * **rpi_params["iconds"]** ( list of ints ): The list of initial geometry indices of the trajectory segment used for setting model_params. 
                Each initial condition in this list characterizes each batch.

            * **rpi_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_params["npatches"]** ( int ): The number of patches. The time duration of each patch dynamics will be int(nsteps/npatches) + 1. 
                The additional single step is because the tsh recipe (see libra_py.dynamics.tsh.compute for details) used for the patch dynamics 
                saves the dynamics information before the electron-nuclear propagation in each MD loop.

            * **rpi_params["nstates"]** ( int ): The number of electronic states.

            * **rpi_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_params["prefix"]** ( string ): The prefix of directories having the output patch dynamics. 
                The patch dynamics data are characterized by three numbers - the batch index, the patch index, the initial state index.
                Thus, the full name of output directory is expressed by these three numbers and their limits as follows.
                  rpi_params["prefix"] + F"n{len(rpi_params["iconds"])}" + "_ibatch{X}" + _n{rpi_params["npatches"]} + F"_ipatch{Y}" \
                        + F"_n{rpi_params["nstates"]}" + F"_istate{Z}"
                  X = 0 to len(rpi_params["iconds"]) - 1; Y = 0 to rpi_params["npatches"] - 1; Z = 0 to rpi_params["nstates"] - 1
                An NBRA trajectory from an initial geometry defined in rpi_params['iconds'] is called a batch. 
        
        The rest of parameters are for setting patch dynamics conducted by `libra_py.dynamics.tsh.compute.generic_recipe`.

        _dyn_params ( dictionary ): control parameters for running the dynamics
            see the documentation of the :func:`libra_py.dynamics.tsh.run_dynamics` function

        compute_model ( PyObject ): the pointer to the Python function that performs the Hamiltonian calculations

        _model_params ( dictionary ): contains the selection of a model and the parameters
            for that model Hamiltonian. In addition to all the parameters, should contain the key

            * **_model_params["model0"]** ( int ): the selection of the model to be used to compute diabatic-to-adiabatic
                transformation. The particular value chosen here depends on the way the `compute_model` functions is defined

               :Note: the function selected should be able to generate the transformation matrix.

        _init_elec ( dictionary ): control parameters for initialization of electronic variables
            see the documentation of the :func:`libra_py.dynamics.init_electronic_dyn_var` function

        _init_nucl (dict): controls how to sample nuclear DOFs from what the user provides. Can contain:
          * **_init_nucl["init_type"]** (int) : the type of the sampling:
            - 0 : Coords and momenta are set exactly to the given value
            - 1 : Coords are set, momenta are sampled
            - 2 : Coords are sampled, momenta are set
            - 3 : Both coords and momenta are samples [ default ]
          * **_init_nucl["ndof"]** (int) : the number of nuclear DOFs [default: 1]
          * **_init_nucl["q"]** ( list of `ndof` doubles ): average coordinate of all trajectories, for each nuclear DOF,
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["p"]** ( list of `ndof` doubles ): average momentum of all trajectories, for each nuclear DOF,
                  a.u. of lenth [ default: [0.0] ]
          * **_init_nucl["mass"]** (list of `ndof` doubles): the mass for each nuclear DOF, in a.u. of mass [ default: [2000.0] ]
          * **_init_nucl["force_constant"]** (list of `ndof` doubles) : the force constant of the harmonic oscillator in each
                  nuclear DOF, determins the width of the coordinate and momentum samplings [ default: [0.001] ]

    Return:
        None: but performs the action
    """
   
    rpi_params = dict(_rpi_params)
    model_params = dict(_model_params)
    dyn_params = dict(_dyn_params)
    init_elec = dict(_init_elec)
    init_nucl = dict(_init_nucl)
    
    # ========== RPI params ==========
    critical_params = []
    default_params = {"run_slurm": False, "submit_template": 'submit_template.slm', "submission_exe": 'sbatch',
                      "run_python_file": 'run_template.py', "python_exe": 'python', "nsteps": 1, "npatches": 1, "nstates":2, "iconds": [0], 
                      "dt": 1.0*units.fs2au, "prefix": 'out_'}
    
    comn.check_input(rpi_params, default_params, critical_params)
    
    # ========== Model params ==========
    comn.check_input(model_params, {}, ["model0"])

    # ========== Nuclear params ==========
    # First check the inputs for the initialization of nuclear variables
    comn.check_input(init_nucl, {"init_type": 3, "ndof": 1, "q": [0.0], "p": [0.0], "mass": [2000.0],
                                 "force_constant": [0.01], "q_width": [1.0], "p_width": [1.0]}, [])
    ndof = len(init_nucl["mass"])

    # ========== Dynamics params ==========
    # Now, we can use the information of ndof to define the default quantum_dofs
    comn.check_input(dyn_params, {"rep_tdse": 1, "rep_sh": 1, "is_nbra": 0, "direct_init": 0, "ntraj": 1,
                                  "quantum_dofs": list(range(ndof))
                                  }, [])

    # RPI setting
    nsteps, dt = rpi_params["nsteps"], rpi_params["dt"]
    iconds, npatches, nstates = rpi_params["iconds"], rpi_params["npatches"], rpi_params["nstates"]
    prefix = rpi_params["prefix"]

    nsteps_patch = int(nsteps / npatches) + 1

    for ibatch, icond in enumerate(iconds):
        for ipatch in range(npatches):
            for istate in range(nstates):
                dir_patch = prefix + F"n{len(iconds)}_ibatch{ibatch}_n{npatches}_ipatch{ipatch}_n{nstates}_istate{istate}"
                if not os.path.exists(dir_patch):
                    os.mkdir(dir_patch)
                os.chdir(dir_patch)

                # Submit jobs or run each patch dynamics through this loop        
                if rpi_params["run_slurm"]:
                    print(F"Submitting a patch dynamics job of ibatch = {ibatch}, ipatch = {ipatch}, istate = {istate}")
                    os.system('cp ../%s %s' % (rpi_params["run_python_file"], rpi_params["run_python_file"]))
                    os.system('cp ../%s %s' % (rpi_params["submit_template"], rpi_params["submit_template"]))
                    
                    icond1 = icond + ipatch * int(nsteps / npatches)
                    
                    file = open(rpi_params["submit_template"], 'a')
                    args_fmt = ' --nsteps %d --icond %d --istate %d --dt %f'
                    file.write('%s %s' % (rpi_params["python_exe"], rpi_params["run_python_file"]) + args_fmt % (nsteps_patch, icond1, istate, dt))
                    file.close()

                    os.system('%s %s' % (rpi_params["submission_exe"], rpi_params["submit_template"]))
                else:
                    print(F"Running patch dynamics of ibatch = {ibatch}, ipatch = {ipatch}, istate = {istate}")
                    run_patch_dyn_serial(rpi_params, dyn_params, compute_model, model_params, init_elec, init_nucl, ibatch, ipatch, istate)
                os.chdir('../')

def print_pop(outfile, time, pops):
    """
    This function prints the RPI population
    """
    with open(outfile, "w") as f:
        for istep in range(time.shape[0]):
            line = f"{time[istep]:13.8f}  " + "  ".join(f"{pops[istep, ist]:13.8f}" for ist in range(pops.shape[1])) + "\n"
            f.write(line)

def process_batch(args):
    """
    This function is used for a parallel execution of the patch summation in the run_sum_rpi function. 
    A single run of this function gives the RPI population on one batch.

    Args:
        input `args` list contains the following:

        ibatch ( int ): the batch index

        rpi_sum_params ( dictionary ): parameters controlling the patch summation in the RPI calculation
            Can contain:
    
            * **rpi_sum_params["nprocs"]** ( int ): The number of processors to be used.

            * **rpi_sum_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_sum_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_sum_params["istate"]** ( int ): The initial state

            * **rpi_sum_params["nbatches"]** ( int ): The number of batches, i.e., the number of initial geometries you used in the previous patch dynamics calculation.
            This corresponds to `len(rpi_params["iconds"])`, where `rpi_params` is a previous patch-dynamics params.

            * **rpi_sum_params["npatches"]** ( int ): The number of patches.

            * **rpi_sum_params["nstates"]** ( int ): The number of electronic states.
            
            * **rpi_sum_params["path_to_save_patch"]** ( string ): The path of the precomputed patch dynamics
    
            * **rpi_sum_params["prefix_patch"]** ( string ): The prefix of patch dynamics directories

            * **rpi_sum_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        pops (nparray): a single-batch population numpy array
    """

    ibatch, rpi_sum_params = args
    
    critical_params = ["path_to_save_patch"]
    default_params = {"nprocs": 1, "nsteps": 1, "dt": 1.0*units.fs2au, "istate": 0, "nbatches": 1, "npatches": 1, "nstates": 2, 
                      "prefix_patch": 'out_', "prefix": 'POP_NPATCH_'}

    comn.check_input(rpi_sum_params, default_params, critical_params)
    
    nsteps, dt, istate = rpi_sum_params["nsteps"], rpi_sum_params["dt"], rpi_sum_params["istate"]
    nbatches, npatches, nstates = rpi_sum_params["nbatches"], rpi_sum_params["npatches"], rpi_sum_params["nstates"]
    path_to_save_patch, prefix_patch = rpi_sum_params["path_to_save_patch"], rpi_sum_params["prefix_patch"]
    prefix = rpi_sum_params["prefix"]
    
    nstep_patch = int(nsteps / npatches)
    pops = np.zeros((npatches * nstep_patch + 1, nstates))
    
    P_temp = np.zeros(nstates)
    P_temp[istate] = 1.0
    pops[0, :] = P_temp
    
    for ipatch in range(npatches):
        print(F"Summing ipatch = {ipatch}, ibatch = {ibatch}")
       
        # Compute the transition probability from a patch dynamics 
        T = np.zeros((nstep_patch, nstates, nstates)) # (timestep in a patch, init, dest)
        for ist in range(nstates):
            with h5py.File(F"{path_to_save_patch}/{prefix_patch}" +
                           F"n{nbatches}_ibatch{ibatch}_n{npatches}_ipatch{ipatch}_n{nstates}_istate{ist}/mem_data.hdf", 'r') as f:
                pop_adi_data = np.array(f["se_pop_adi/data"])
            T[:, ist, :] = pop_adi_data[1:,:]
            T[:, ist, ist] = 0.0  # zero diagonal explicitly 
       
        # Fill diagonal by the norm conservation
        T[:, np.arange(nstates), np.arange(nstates)] = 1.0 - np.sum(T, axis=2)
        
        for j in range(nstep_patch):
            glob_time = ipatch * nstep_patch + j + 1
            pops[glob_time, :] = P_temp @ T[j]
            if j == nstep_patch - 1:
                P_temp = pops[glob_time]
    
    # Per-batch output
    time = np.array([x for x in range(npatches * nstep_patch + 1)]) * dt
    print(F"Print the population from ibatch")
    print_pop(prefix + F"{npatches}_ibatch{ibatch}.dat", time, pops)
    
    return pops

def run_sum(rpi_sum_params):
    """
    This function conducts the RPI patch summation to yield the population dynamics in the whole time domain.

    Args:

        rpi_sum_params ( dictionary ): parameters controlling the patch summation in the RPI calculation
            Can contain:
    
            * **rpi_sum_params["nprocs"]** ( int ): The number of processors to be used.

            * **rpi_sum_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_sum_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_sum_params["istate"]** ( int ): The initial state

            * **rpi_sum_params["nbatches"]** ( int ): The number of batches, i.e., the number of initial geometries you used in the previous patch dynamics calculation.
            This corresponds to `len(rpi_params["iconds"])`, where `rpi_params` is a previous patch-dynamics params.

            * **rpi_sum_params["npatches"]** ( int ): The number of patches.

            * **rpi_sum_params["nstates"]** ( int ): The number of electronic states.
            
            * **rpi_sum_params["path_to_save_patch"]** ( string ): The path of the precomputed patch dynamics
    
            * **rpi_sum_params["prefix_patch"]** ( string ): The prefix of patch dynamics directories
            
            * **rpi_sum_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        None: but performs the action
    """
    critical_params = ["path_to_save_patch"]
    default_params = {"nprocs": 1, "nsteps": 1, "dt": 1.0*units.fs2au, "istate": 0, "nbatches": 1, "npatches": 1, "nstates": 2, 
                      "prefix_patch": 'out_', "prefix": 'POP_NPATCH_'}

    comn.check_input(rpi_sum_params, default_params, critical_params)

    nprocs = rpi_sum_params["nprocs"]
  
    nsteps, dt = rpi_sum_params["nsteps"], rpi_sum_params["dt"]
    nbatches, npatches, nstates = rpi_sum_params["nbatches"], rpi_sum_params["npatches"], rpi_sum_params["nstates"]
    prefix = rpi_sum_params["prefix"]

    nstep_patch = int(nsteps / npatches)
    time = np.array([x for x in range(npatches * nstep_patch + 1)]) * dt
    pops_avg = np.zeros((npatches * nstep_patch + 1, nstates))
    
    with mp.Pool(processes=nprocs) as pool:
        results = pool.map(process_batch, [(ibatch, rpi_sum_params) for ibatch in range(nbatches)])
    
    for pops in results:
        pops_avg += pops
    
    pops_avg /= nbatches
    
    print("Print the final population from all batches")
    print_pop(prefix + F"{npatches}_all.dat", time, pops_avg)


def run_sum_crude(rpi_sum_params):
    """
    This function conducts the RPI patch summation to yield the population dynamics in the whole time domain.

    Args:
        
        rpi_sum_params ( dictionary ): parameters controlling the patch summation in the RPI calculation
            Can contain:
    
            * **rpi_sum_params["nprocs"]** ( int ): The number of processors to be used.

            * **rpi_sum_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_sum_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_sum_params["istate"]** ( int ): The initial state

            * **rpi_sum_params["nbatches"]** ( int ): The number of batches, i.e., the number of initial geometries you used in the previous patch dynamics calculation.
            This corresponds to `len(rpi_params["iconds"])`, where `rpi_params` is a previous patch-dynamics params.

            * **rpi_sum_params["npatches"]** ( int ): The number of patches.

            * **rpi_sum_params["nstates"]** ( int ): The number of electronic states.
            
            * **rpi_sum_params["path_to_save_patch"]** ( string ): The path of the precomputed patch dynamics
    
            * **rpi_sum_params["prefix_patch"]** ( string ): The prefix of patch dynamics directories
            
            * **rpi_sum_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        None: but performs the action
    """
    critical_params = ["path_to_save_patch"]
    default_params = {"nprocs": 1, "nsteps": 1, "dt": 1.0*units.fs2au, "istate": 0, "nbatches": 1, "npatches": 1, "nstates": 2, 
                      "prefix_patch": 'out_', "prefix": 'POP_NPATCH_'}

    comn.check_input(rpi_sum_params, default_params, critical_params)

    nsteps, dt, istate = rpi_sum_params["nsteps"], rpi_sum_params["dt"], rpi_sum_params["istate"]

    nbatches, npatches, nstates = rpi_sum_params["nbatches"], rpi_sum_params["npatches"], rpi_sum_params["nstates"]
    path_to_save_patch, prefix_patch = rpi_sum_params["path_to_save_patch"], rpi_sum_params["prefix_patch"]
    prefix = rpi_sum_params["prefix"]

    nstep_patch = int(nsteps / npatches)

    time = np.array([x for x in range(npatches * nstep_patch + 1)]) * dt

    # The total population across all batches
    pops_avg = np.zeros((npatches * nstep_patch + 1, nstates))
    
    for ibatch in range(nbatches):
        # The global population on a batch
        pops = np.zeros((npatches * nstep_patch + 1, nstates))

        P_temp = np.zeros(nstates)
        P_temp[istate] = 1.0  # Set the initial population
        pops[0, :] = P_temp
   
        pop_patch = np.zeros((nstep_patch, nstates, nstates)) # A temporary array to read population of a patch dynamics
        T = np.zeros((nstep_patch, nstates, nstates)) # The transition probability from a patch dynamics
        for ipatch in range(npatches):
            print(F"Summing ipatch = {ipatch}, ibatch = {ibatch}")
            pop_patch.fill(0.0)
            T.fill(0.0)
            for ist in range(nstates):
                with h5py.File(F"{path_to_save_patch}/{prefix_patch}" +
                               F"n{nbatches}_ibatch{ibatch}_n{npatches}_ipatch{ipatch}_n{nstates}_istate{ist}/mem_data.hdf", 'r') as f:
                    pop_adi_data = np.array(f["se_pop_adi/data"])
            
                for istep in range(nstep_patch):
                    for jst in range(nstates):
                        pop_patch[istep, jst, ist] = pop_adi_data[1 + istep, jst]

            for ist in range(nstates):
                for jst in range(nstates):
                    if ist == jst:
                        continue
                    T[:, ist, jst] = pop_patch[:, jst, ist]

            for istep in range(nstep_patch):
                for ist in range(nstates):
                    T[istep, ist, ist] = 1. - np.sum(T[istep, ist, :])

            for j in range(nstep_patch):
                glob_time = ipatch * nstep_patch + j + 1
                for jst in range(nstates):
                    pops[glob_time, jst] = np.sum(P_temp * T[j, :, jst])
                if j == nstep_patch - 1:
                    P_temp = pops[glob_time]
        pops_avg += pops
    
        print(F"Print the population from ibatch, ibatch = {ibatch}")
        print_pop(prefix + F"{npatches}_ibatch{ibatch}.dat", time, pops)

    pops_avg /= nbatches

    print("Print the final population from all batches")
    print_pop(prefix + F"{npatches}_all.dat", time, pops_avg)

