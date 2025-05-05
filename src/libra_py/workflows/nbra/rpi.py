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
    
class tmp:
    pass

def compute_model(q, params, full_id):
    """
    This function serves as an interface function for a serial patch dynamics calculation.
    """

    timestep = params["timestep"]
    nst = params["nstates"]
    E = params["E"]
    NAC = params["NAC"]
    Hvib = params["Hvib"]
    St = params["St"]

    obj = tmp()

    obj.ham_adi = data_conv.nparray2CMATRIX( np.diag(E[timestep, : ]) )
    obj.nac_adi = data_conv.nparray2CMATRIX( NAC[timestep, :, :] )
    obj.hvib_adi = data_conv.nparray2CMATRIX( Hvib[timestep, :, :] )
    obj.basis_transform = CMATRIX(nst,nst); obj.basis_transform.identity()  #basis_transform
    obj.time_overlap_adi = data_conv.nparray2CMATRIX( St[timestep, :, :] )

    return obj


def run_patch_dyn_serial(rpi_params, ibatch, ipatch, istate):
    """
    This function runs a single patch dynamics calculation for the time segment from istep to fstep and a given initial state.
    For conducting a coherent Ehrenfest propagation, the `libra_py.dynamics.tsh.compute` module is used.

    Args:

        rpi_params ( dictionary ): parameters controlling the patch dynamics in the RPI calculation
            Can contain:

            * **rpi_params["run_slurm"]** ( bool ): Whether to use the slurm environment to submit the jobs using the submit_template file. 
                If it is set to False, it will run the calculations on the active session but multiple jobs will be run on the current active session.

            * **rpi_params["submit_template"]** ( string ): The path of a template slurm submit file.

            * **rpi_params["submission_exe"]** ( string ): The submission executable

            * **rpi_params["run_python_file"]** ( string ): The path of a template running script.
            
            * **rpi_params["python_exe"]** ( string ): The Python executable

            * **rpi_params["path_to_save_Hvibs"]** ( string ): The path of the vibronic Hamiltonian files.
            
            * **rpi_params["basis_type"]** ( string ): The electronic basis type (such as the Slater determinant or CI adpatation, i.e., 'sd' or 'ci').

            * **rpi_params["iread"]** ( int ): The initial step to read the vibronic Hamiltonian and time overlap

            * **rpi_params["fread"]** ( int ): The final step to read the vibronic Hamiltonian and time overlap

            * **rpi_params["iconds"]** ( list of ints ): The list of initial step indices from the trajectory segment from iread to fread. 
                Each initial condition characterizes each batch.

            * **rpi_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_params["npatches"]** ( int ): The number of patches. The time duration of each patch dynamics will be int( nsteps / npatches ) + 1. 
              The additional single step is because the tsh recipe (see libra_py.dynamics.tsh.compute for details) used for the patch dynamics 
              saves the dynamics information before the electron-nuclear propagation in each MD loop.

            * **rpi_params["nstates"]** ( int ): The number of electronic states.

            * **rpi_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_params["path_to_save_patch"]** ( string ): The path of the output patch dynamics

        ibatch ( int ): the initial geometry index of rpi_params["iconds"]
        
        ipatch ( int ): the patch index, running from 0 to rpi_params["npatches"] - 1

        istate ( int ): the initial electronic state for this patch dynamics

    Return:
        None: but performs the action
    """
    
    critical_params = ["path_to_save_patch"]
    default_params = {"run_slurm": False, "submit_template": 'submit_template.slm', "submission_exe": 'sbatch',
                      "run_python_file": 'run_template.py', "python_exe": 'python', "path_to_save_Hvibs": 'res', "basis_type": 'ci', 
                      "iread": 0, "fread": 1, "nsteps": 1, "npatches": 1, "nstates":2, "iconds": [0], "dt": 1.0*units.fs2au}
    
    comn.check_input(rpi_params, default_params, critical_params)
    
    path_to_save_Hvibs, basis_type = rpi_params["path_to_save_Hvibs"], rpi_params["basis_type"]
    iread, fread = rpi_params["iread"], rpi_params["fread"]
    nsteps, dt = rpi_params["nsteps"], rpi_params["dt"]
    iconds, npatches, nstates = rpi_params["iconds"], rpi_params["npatches"], rpi_params["nstates"]

    # Compute the istep and fstep here, for each patch
    istep = iread + iconds[ibatch] + ipatch * int( nsteps / npatches )
    fstep = istep + int( nsteps / npatches ) + 1
    
    NSTEPS = fstep - istep # The patch duration

    # Set the NBRA model by reading the vibronic Hamiltonian and time overlap data
    model_params = {"timestep":0, "icond":0,  "model0":0, "nstates":nstates}
    
    E = []
    for step in range(istep, fstep):
        energy_filename = F"{path_to_save_Hvibs}/Hvib_{basis_type}_{step}_re.npz"
        energy_mat = sp.load_npz(energy_filename)
        E.append( np.array( np.diag( energy_mat.todense() ) ) )
    E = np.array(E)
    
    St = []
    for step in range(istep, fstep):
        St_filename = F"{path_to_save_Hvibs}/St_{basis_type}_{step}_re.npz"
        St_mat = sp.load_npz(St_filename)
        St.append( np.array( St_mat.todense() ) )
    St = np.array(St)

    NAC = []
    Hvib = []
    for c, step in enumerate(range(istep, fstep)):
        nac_filename = F"{path_to_save_Hvibs}/Hvib_{basis_type}_{step}_im.npz"
        nac_mat = sp.load_npz(nac_filename)
        NAC.append( np.array( nac_mat.todense() ) )
        Hvib.append( np.diag(E[c, :])*(1.0+1j*0.0)  - (0.0+1j)*nac_mat[:, :] )
    
    NAC = np.array(NAC)
    Hvib = np.array(Hvib)

    model_params.update({"E": E, "St": St, "NAC": NAC, "Hvib": Hvib})
    
    # Setting the coherent Ehrenfest propagation. Define the argument-dependent part first.
    dyn_general = {"nsteps":NSTEPS, "nstates":nstates, "dt":dt, "nfiles": NSTEPS, "which_adi_states":range(nstates), "which_dia_states":range(nstates)}
    
    dyn_general.update({"ntraj":1, "mem_output_level":2, "progress_frequency":1, "properties_to_save":["timestep", "time", "se_pop_adi"], "prefix":"out", "isNBRA":0, 
                   "ham_update_method":2, "ham_transform_method":0, "time_overlap_method":0, "nac_update_method":0,
                   "hvib_update_method":0, "force_method":0, "rep_force":1, "hop_acceptance_algo":0, "momenta_rescaling_algo":0, "rep_tdse":1,
                   "electronic_integrator":2, "tsh_method":-1, "decoherence_algo":-1, "decoherence_times_type":-1, "do_ssy":0, "dephasing_informed":0})

    # Nuclear DOF - these parameters don't matter much in the NBRA calculations
    nucl_params = {"ndof":1, "init_type":3, "q":[-10.0], "p":[0.0], "mass":[2000.0], "force_constant":[0.0], "verbosity":-1 }

    # Electronic DOF
    elec_params = {"ndia":nstates, "nadi":nstates, "verbosity":-1, "init_dm_type":0, "init_type":1, "rep":1,"istate":istate }

    rnd = Random()

    res = tsh_dynamics.generic_recipe(dyn_general, compute_model, model_params, elec_params, nucl_params, rnd)

def run_patch_rpi(rpi_params):
    """
    This function distributes the jobs to perform patch dynamics calculations in the restricted path integral (RPI) method.

    Args:

        rpi_params ( dictionary ): parameters controlling the patch dynamics in the RPI calculation
            Can contain:

            * **rpi_params["run_slurm"]** ( bool ): Whether to use the slurm environment to submit the jobs using the submit_template file. 
                If it is set to False, it will run the calculations on the active session but multiple jobs will be run on the current active session.

            * **rpi_params["submit_template"]** ( string ): The path of a template slurm submit file.

            * **rpi_params["submission_exe"]** ( string ): The submission executable

            * **rpi_params["run_python_file"]** ( string ): The path of a template running script.
            
            * **rpi_params["python_exe"]** ( string ): The Python executable

            * **rpi_params["path_to_save_Hvibs"]** ( string ): The path of the vibronic Hamiltonian files.
            
            * **rpi_params["basis_type"]** ( string ): The electronic basis type (such as the Slater determinant or CI adpatation, i.e., 'sd' or 'ci').

            * **rpi_params["iread"]** ( int ): The initial step to read the vibronic Hamiltonian and time overlap

            * **rpi_params["fread"]** ( int ): The final step to read the vibronic Hamiltonian and time overlap

            * **rpi_params["iconds"]** ( list of ints ): The list of initial step indices from the trajectory segment from iread to fread. 
                Each initial condition characterizes each batch.

            * **rpi_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_params["npatches"]** ( int ): The number of patches. The time duration of each patch dynamics will be int(nsteps/npatches) + 1. 
              The additional single step is because the tsh recipe (see libra_py.dynamics.tsh.compute for details) used for the patch dynamics 
              saves the dynamics information before the electron-nuclear propagation in each MD loop.

            * **rpi_params["nstates"]** ( int ): The number of electronic states.

            * **rpi_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_params["path_to_save_patch"]** ( string ): The path of the output patch dynamics

    Return:
        None: but performs the action
    """
    
    critical_params = ["path_to_save_patch"]
    default_params = {"run_slurm": False, "submit_template": 'submit_template.slm', "submission_exe": 'sbatch',
                      "run_python_file": 'run_template.py', "python_exe": 'python', "path_to_save_Hvibs": 'res', "basis_type": 'ci', 
                      "iread": 0, "fread": 1, "nsteps": 1, "npatches": 1, "nstates":2, "iconds": [0], "dt": 1.0*units.fs2au}
    
    comn.check_input(rpi_params, default_params, critical_params)

    path_to_save_Hvibs, basis_type = rpi_params["path_to_save_Hvibs"], rpi_params["basis_type"]
    iread, fread = rpi_params["iread"], rpi_params["fread"]
    nsteps, dt = rpi_params["nsteps"], rpi_params["dt"]
    iconds, npatches, nstates = rpi_params["iconds"], rpi_params["npatches"], rpi_params["nstates"]

    out_dir = rpi_params["path_to_save_patch"]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    os.chdir(out_dir)

    for ibatch, icond in enumerate(iconds):
        for ipatch in range(npatches):
            for istate in range(nstates):
                
                # Compute the istep and fstep here, for each patch
                istep = iread + icond + ipatch * int( nsteps / npatches )
                fstep = istep + int( nsteps / npatches ) + 1

                dir_patch = F"job_{ibatch}_{ipatch}_{istate}"
                if os.path.exists(dir_patch):
                    os.system('rm -rf ' + dir_patch)
                os.mkdir(dir_patch)
                os.chdir(dir_patch)

                # Submit jobs or run each patch dynamics through this loop        
                if rpi_params["run_slurm"]:
                    print(F"Submitting a patch dynamics job of icond = {ibatch}, ipatch = {ipatch}, istate = {istate}")
                    os.system('cp ../../%s %s' % (rpi_params["run_python_file"], rpi_params["run_python_file"]))
                    os.system('cp ../../%s %s' % (rpi_params["submit_template"], rpi_params["submit_template"]))
                    
                    file = open(rpi_params["submit_template"], 'a')
                    args_fmt = ' --path_to_save_Hvibs %s --basis_type %s --istep %d --fstep %d --dt %f --istate %d'
                    file.write('%s %s' % (rpi_params["python_exe"], rpi_params["run_python_file"]) + \
                                args_fmt % (path_to_save_Hvibs, basis_type, istep, fstep, dt, istate) + " >log" )
                    file.close()

                    os.system('%s %s' % (rpi_params["submission_exe"], rpi_params["submit_template"]))
                else:
                    print(F"Running patch dynamics of icond = {ibatch}, ipatch = {ipatch}, istate = {istate}")
                    run_patch_dyn_serial(rpi_params, ibatch, ipatch, istate)
                os.chdir('../')

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
            
            * **rpi_sum_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        pops (nparray): a single-batch population numpy array
    """

    ibatch, rpi_sum_params = args
    
    critical_params = ["path_to_save_patch"]
    default_params = {"nprocs": 1, "nsteps": 1, "dt": 1.0*units.fs2au, "istate": 0, "nbatches": 1, "npatches": 1, "nstates": 2, "prefix": 'out'}

    comn.check_input(rpi_sum_params, default_params, critical_params)
    
    nsteps, dt, istate = rpi_sum_params["nsteps"], rpi_sum_params["dt"], rpi_sum_params["istate"]
    npatches, nstates = rpi_sum_params["npatches"], rpi_sum_params["nstates"]
    path_to_save_patch = rpi_sum_params["path_to_save_patch"]
    
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
            with h5py.File(F"{path_to_save_patch}/job_{ibatch}_{ipatch}_{ist}/out/mem_data.hdf", 'r') as f:
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
    print_pop(rpi_params["prefix"] + F"_ibatch{ibatch}.dat", time, pops)
    
    return pops

def run_sum_rpi(rpi_sum_params):
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
            
            * **rpi_sum_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        None: but performs the action
    """
    critical_params = ["path_to_save_patch"]
    default_params = {"nprocs": 1, "nsteps": 1, "dt": 1.0*units.fs2au, "istate": 0, "nbatches": 1, "npatches": 1, "nstates": 2, "prefix": 'out'}

    comn.check_input(rpi_sum_params, default_params, critical_params)

    nprocs = rpi_sum_params["nprocs"]
  
    nsteps, dt = rpi_sum_params["nsteps"], rpi_sum_params["dt"]
    nbatches, npatches, nstates = rpi_sum_params["nbatches"], rpi_sum_params["npatches"], rpi_sum_params["nstates"]

    nstep_patch = int(nsteps / npatches)
    time = np.array([x for x in range(npatches * nstep_patch + 1)]) * dt
    pops_avg = np.zeros((npatches * nstep_patch + 1, nstates))
    
    with mp.Pool(processes=nprocs) as pool:
        results = pool.map(process_batch, [(ibatch, rpi_sum_params) for ibatch in range(nbatches)])
    
    for pops in results:
        pops_avg += pops
    
    pops_avg /= nbatches
    
    print("Print the final population from all batches")
    print_pop(rpi_sum_params["prefix"] + "_all.dat", time, pops_avg)


def run_sum_rpi_crude(rpi_sum_params):
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
            
            * **rpi_sum_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        None: but performs the action
    """
    critical_params = ["path_to_save_patch"]
    default_params = {"nprocs": 1, "nsteps": 1, "dt": 1.0*units.fs2au, "istate": 0, "nbatches": 1, "npatches": 1, "nstates": 2, "prefix": 'out'}

    comn.check_input(rpi_sum_params, default_params, critical_params)

    nsteps, dt, istate = rpi_sum_params["nsteps"], rpi_sum_params["dt"], rpi_sum_params["istate"]

    nbatches, npatches, nstates = rpi_sum_params["nbatches"], rpi_sum_params["npatches"], rpi_sum_params["nstates"]
    path_to_save_patch = rpi_sum_params["path_to_save_patch"]

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
                with h5py.File(F"{path_to_save_patch}/job_{ibatch}_{ipatch}_{ist}/out/mem_data.hdf", 'r') as f:
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
        print_pop(rpi_sum_params["prefix"] + F"_ibatch{ibatch}.dat", time, pops)

    pops_avg /= nbatches

    print("Print the final population from all batches")
    print_pop(rpi_sum_params["prefix"] + "_all.dat", time, pops_avg)

