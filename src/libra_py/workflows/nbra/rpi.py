import libra_py.units as units
from libra_py import data_conv
import libra_py.data_read as data_read
import util.libutil as comn
import matplotlib.pyplot as plt
import sys
import cmath
import math
import os
import multiprocessing as mp
import time
import numpy as np
import h5py
import scipy.sparse as sp
import multiprocessing as mp

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

def run_patch_rpi(rpi_params):
    """
    This function distributes the jobs to perform patch dynamics calculations in the restricted path integral (RPI) method.

    Args:

        rpi_params ( dictionary ): parameters controlling the execution of the RPI dynamics
            Can contain:

            * **rpi_params["run_slurm"]** ( bool ): Whether to use the slurm environment to submit the jobs using the submit_template file. 
                If it is set to False, it will run the calculations on the active session but multiple jobs will be run on the current active session.

            * **rpi_params["submit_template"]** ( string ): The path of a template slurm submit file.

            * **rpi_params["run_python_file"]** ( string ): The path of a template running script.

            * **rpi_params["path_to_save_Hvibs"]** ( string ): The path of the vibronic Hamiltonian files.

            * **rpi_params["submission_exe"]** ( string ): The submission executable

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

    out_dir = rpi_params["path_to_save_patch"]
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    file = open(rpi_params["run_python_file"], 'r')
    lines = file.readlines()
    file.close()

    os.chdir(out_dir)

    for ibatch, icond in enumerate(rpi_params["iconds"]):
        for ipatch in range(rpi_params["npatches"]):
            for istate in range(rpi_params["nstates"]):
                print(F"Submitting patch dynamics job of icond = {ibatch}, ipatch = {ipatch}, istate = {istate}")
                
                # Compute the istep and fstep here, for each patch
                istep = rpi_params['iread'] + icond + ipatch*int(rpi_params['nsteps']/rpi_params['npatches'])
                fstep = fstep = istep + int(rpi_params['nsteps']/rpi_params['npatches']) + 1

                dir_patch = F"job_{ibatch}_{ipatch}_{istate}"
                if os.path.exists(dir_patch):
                    os.system('rm -rf ' + dir_patch)
                os.mkdir(dir_patch)
                os.chdir(dir_patch)
                file = open('run.py', 'w')

                for i in range(len(lines)):
                    if "params['path_to_save_Hvibs'] =" in lines[i]:
                        file.write("""params['path_to_save_Hvibs'] = "%s"\n""" % rpi_params["path_to_save_Hvibs"])
                    elif "params['istep'] =" in lines[i]:
                        file.write("params['istep'] = %d\n" % istep)
                    elif "params['fstep'] =" in lines[i]:
                        file.write("params['fstep'] = %d\n" % fstep)
                    elif "params['nsteps'] =" in lines[i]:
                        file.write("params['nsteps'] = %d\n" % rpi_params["nsteps"])
                    elif "params['npatches'] =" in lines[i]:
                        file.write("params['npatches'] = %d\n" % rpi_params["npatches"])
                    elif "params['istate'] =" in lines[i]:
                        file.write("params['istate'] = %d\n" % istate)
                    elif "params['dt'] =" in lines[i]:
                        file.write("params['dt'] = %f\n" % rpi_params["dt"])
                    else:
                        file.write(lines[i])
                file.close()
        
                if rpi_params["run_slurm"]:
                    os.system('cp ../../%s %s' % (rpi_params["submit_template"], rpi_params["submit_template"]))
                    os.system('%s %s' % (rpi_params["submission_exe"], rpi_params["submit_template"]))
                else:
                    # Just in case you want to use a bash file and not submitting
                    os.system("python run.py > log")
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
    ibatch, icond, rpi_params = args
    nsteps, dt, istate = rpi_params["nsteps"], rpi_params["dt"], rpi_params["istate"]
    npatches, nstates = rpi_params["npatches"], rpi_params["nstates"]
    path_to_save_patch = rpi_params["path_to_save_patch"]
    
    nstep_patch = int(nsteps / npatches)
    pops = np.zeros((npatches * nstep_patch + 1, nstates))
    
    P_temp = np.zeros(nstates)
    P_temp[istate] = 1.0
    pops[0, :] = P_temp
    
    for ipatch in range(npatches):
        print(F"Summing ipatch = {ipatch}, ibatch = {ibatch}, icond = {icond}")
       
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
    print(F"Print the population from ibatch, icond = {ibatch}, {icond}")
    print_pop(rpi_params["prefix"] + F"_ibatch{ibatch}.dat", time, pops)
    
    return pops

def run_sum_rpi(rpi_params):
    """
    This function conducts the RPI patch summation to yield the population dynamics in the whole time domain.

    Args:

        rpi_params ( dictionary ): parameters controlling the execution of the RPI dynamics
            Can contain:
    
            * **rpi_params["nprocs"]** ( int ): The number of processors to be used.

            * **rpi_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_params["istate"]** ( int ): The initial state

            * **rpi_params["iconds"]** ( list of ints ): The list of initial step indices from the trajectory segment from iread to fread. 
                Each initial condition characterizes each batch.

            * **rpi_params["npatches"]** ( int ): The number of patches.

            * **rpi_params["nstates"]** ( int ): The number of electronic states.
            
            * **rpi_params["path_to_save_patch"]** ( string ): The path of the precomputed patch dynamics
            
            * **rpi_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        None: but performs the action
    """
    nprocs = rpi_params["nprocs"]
  
    nsteps, dt = rpi_params["nsteps"], rpi_params["dt"]
    npatches, nstates = rpi_params["npatches"], rpi_params["nstates"]
    iconds = rpi_params["iconds"]

    nstep_patch = int(nsteps / npatches)
    time = np.array([x for x in range(npatches * nstep_patch + 1)]) * dt
    pops_avg = np.zeros((npatches * nstep_patch + 1, nstates))
    
    with mp.Pool(processes=nprocs) as pool:
        results = pool.map(process_batch, [(ibatch, icond, rpi_params) for ibatch, icond in enumerate(iconds)])
    
    for pops in results:
        pops_avg += pops
    
    pops_avg /= len(iconds)
    
    print("Print the final population from all batches")
    print_pop(rpi_params["prefix"] + "_all.dat", time, pops_avg)


def run_sum_rpi_crude(rpi_params):
    """
    This function conducts the RPI patch summation to yield the population dynamics in the whole time domain.

    Args:

        rpi_params ( dictionary ): parameters controlling the execution of the RPI dynamics
            Can contain:
            
            * **rpi_params["nsteps"]** ( int ): The total number of RPI simulation steps

            * **rpi_params["dt"]** ( double ): the time step in the atomic unit.
            
            * **rpi_params["istate"]** ( int ): The initial state

            * **rpi_params["iconds"]** ( list of ints ): The list of initial step indices from the trajectory segment from iread to fread. 
                Each initial condition characterizes each batch.

            * **rpi_params["npatches"]** ( int ): The number of patches.

            * **rpi_params["nstates"]** ( int ): The number of electronic states.
            
            * **rpi_params["path_to_save_patch"]** ( string ): The path of the precomputed patch dynamics
            
            * **rpi_params["prefix"]** ( string ): The prefix for the population dynamics output

    Return:
        None: but performs the action
    """

    nsteps, dt, istate = rpi_params["nsteps"], rpi_params["dt"], rpi_params["istate"]

    iconds, npatches, nstates = rpi_params["iconds"], rpi_params["npatches"], rpi_params["nstates"]
    path_to_save_patch = rpi_params["path_to_save_patch"]

    nstep_patch = int(nsteps / npatches)

    time = np.array([x for x in range(npatches * nstep_patch + 1)]) * dt

    # The total population across all batches
    pops_avg = np.zeros((npatches * nstep_patch + 1, nstates))
    
    for ibatch, icond in enumerate(iconds):
        # The global population on a batch
        pops = np.zeros((npatches * nstep_patch + 1, nstates))

        P_temp = np.zeros(nstates)
        P_temp[istate] = 1.0  # Set the initial population
        pops[0, :] = P_temp
   
        pop_patch = np.zeros((nstep_patch, nstates, nstates)) # A temporary array to read population of a patch dynamics
        T = np.zeros((nstep_patch, nstates, nstates)) # The transition probability from a patch dynamics
        for ipatch in range(npatches):
            print(F"Summing ipatch = {ipatch}, ibatch = {ibatch}, icond = {icond}")
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
    
        print(F"Print the population from ibatch, icond = {ibatch}, {icond}")
        print_pop(rpi_params["prefix"] + F"_ibatch{ibatch}.dat", time, pops)

    pops_avg /= len(iconds)

    print("Print the final population from all batches")
    print_pop(rpi_params["prefix"] + "_all.dat", time, pops_avg)

