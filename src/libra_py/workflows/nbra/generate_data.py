# ***********************************************************
# * Copyright (C) 2024 Mohammad Shakiba and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
# ***********************************************************/

"""
.. module:: generate_data
   :platform: Unix
   :synopsis: This module implements functions for running quantum chemistry softwares
   such as CP2K and DFTB+ for a precomputed trajectory and generate the atomic orbital
   matrices such as Kohn-Sham Hamiltonian, density matrix, overlap matrix, etc.

.. moduleauthor:: Mohammad Shakiba and Alexey V. Akimov

"""

import os
import json
import time
import numpy as np
import sys
import glob
import libra_py.packages.cp2k.methods as CP2K_methods
import util.libutil as comn


def make_input(prefix, input_template, software, trajectory_xyz_file, step):
    """
    This function makes an input for a step
    based on the input template and the trajectory xyz file
    Args:
        prefix (string): The prefix with which the inputs names are saved
        input_template (string): The full path to input template
        software (string): The software name for which the calculations should be done
        trajectory_xyz_file (string): The full path to the trajectory xyz file
        step (integer): The step^th geometry in the molecular dynamics trajectory
                        for which the calculations should be done
    Returns:
        None
    """
    f = open(input_template, 'r')
    lines_input = f.readlines()
    f.close()

    _, _ = CP2K_methods.read_trajectory_xyz_file(trajectory_xyz_file, step)

    if software.lower() == "cp2k":
        f = open(F"input_{prefix}_{step}.inp", "w")
        for i in range(len(lines_input)):
            if "COORD_FILE_NAME".lower() in lines_input[i].lower():
                f.write(F"     COORD_FILE_NAME  coord-{step}.xyz \n")
            elif "project " in lines_input[i].lower() or "project_name" in lines_input[i].lower():
                f.write(F"  PROJECT {prefix}_{step}\n")
            elif "filename" in lines_input[i].lower() and "!" not in lines_input[i]:
                f.write("     FILENAME libra\n")
            else:
                f.write(lines_input[i])
        f.close()
    elif software.lower() == "dftb+":
        f = open(F"dftb_in_{prefix}_{step}.hsd", "w")
        for i in range(len(lines_input)):
            if "geometry" in lines_input[i].lower():
                geometry_line = i
                break
        for i in range(geometry_line, len(lines_input)):
            if "{" in lines_input[i]:
                start = i
            elif "}" in lines_input[i]:
                end = i
                break
        f.write("Geometry = xyzFormat {\n")
        f.write(F"<<< coord-{step}.xyz \n")
        f.write("}\n")
        for i in range(end + 1):
            f.write(lines_input[i])
        f.close()


def save_data(params, prefix, step, directory, guess):
    """
    This is an auxiliary function to extract and save the output matrices data
    Args:
        params (dict): The detailed explanation of the parameters are given in the distribute_jobs function.
        prefix (string): The prefix with which the data should be saved.
        step (integer): The step^th geometry in the molecular dynamics trajectory.
        directory (string): The full path to which the data will be saved.
        guess (bool): An internal boolean flag showing whether the guess or reference calculations should be saved.
    Returns:
        None
    """
    if guess:
        software = params["guess_software"]
    else:
        software = params["reference_software"]

    if software.lower() == "cp2k":
        try:
            filename = F"{prefix}_{step}-libra-1_0.Log"
            properties, data = CP2K_methods.read_ao_matrices(filename)
            for i in range(len(data)):
                tmp_name = "".join(properties[i].split()).replace("-", "_").lower()
                np.save(F"{directory}/{prefix}_{tmp_name}_{step}.npy", data[i])
            if step == params["istep"]:
                # This is how we can get the sample molden, wfn, and Log files.
                os.system("mkdir ../sample_files")
                os.system(f"mv {prefix}_{step}* ../sample_files")
                os.system("cp params.json ../sample_files")
            if params["remove_raw_outputs"]:  # and step>params["istep"]:
                os.system(F"rm {prefix}_{step}*Log {prefix}_{step}*.wfn* {prefix}_{step}*molden")
        except BaseException:
            raise ("Could not read the data from CP2K .Log files!")

    elif software.lower() == "dftb+":
        try:
            filename = "hamsqrd1.dat"
            if os.path.exists(filename):
                data = np.loadtxt(filename, skiprows=5)
                np.save(F"{prefix}_ham_{step}.npy", data)
            else:
                raise ("Could not find file hamsqrd1.dat...")
        except BaseException:
            raise ("Could not read data from DFTB+ .dat files!")


def run_single_point(params, prefix, step, guess):
    """
    This function runs the software on the Linux environment
    Args:
        params (dict): The detailed explanation of the parameters are given in the distribute_jobs function.
        prefix (string): The prefix showing for which input the data should be run.
        step (integer): THe step^th geometry in the molecular dynamics trajectory.
        guess (bool): An internal boolean flag showing whether the guess or reference calculations should be saved.
    Returns:
        None
    """
    if guess:
        software = params["guess_software"]
        software_exe = params["guess_software_exe"]
        mpi_exe = params["guess_mpi_exe"]
    else:
        software = params["reference_software"]
        software_exe = params["reference_software_exe"]
        mpi_exe = params["reference_mpi_exe"]

    if software.lower() == "cp2k":
        os.system(
            F"{mpi_exe} -np {params['nprocs']} {software_exe} -i input_{prefix}_{step}.inp -o output_{prefix}_{step}.log")
    elif software.lower() == "dftb+":
        os.system(F"export OMP_NUM_THREADS={nprocs}")
        os.system(F"mv dftb_in_{prefix}_{step}.hsd dftb_in.hsd")
        os.system(F"{software_exe} dftb_in.hsd > output_{prefix}_{i}.log")
        os.system(F"mv dftb_in.hsd dftb_in_{prefix}_{step}.hsd")


def gen_data(params, step):
    """
    This function takes multiple parameters to run electronic structure calculations
    Args:
        params (dict): The detailed explanation of the parameters are given in the distribute_jobs function.
        step (integer): The step^th geometry of the molecular dynamics trajectory
    Returns:
        None
    """
    if params["do_guess"]:
        # Addition on 4/23/2024: Saving the guess Hams for a full trajectory is not wise since
        # for some systems this might overflow the disk. Instead, we will recompute these cheap
        # guesses when we want to compute the properties after the model is trained in the ml_map module
        if step in params["reference_steps"]:
            # Generate the step^th geometry from the trajectory file
            # and create an input file for that
            make_input(
                params["prefix"],
                params["guess_input_template"],
                params["guess_software"],
                params["trajectory_xyz_file"],
                step)
            # Run the calculations for this input file
            t1 = time.time()
            print("Guess calculations for step", step)
            run_single_point(params, params["prefix"], step, True)

            save_data(params, params["prefix"], step, params["guess_dir"], True)
            # else:
            #    print("to be implemented!")
            print(f"Elapsed time for guess calculations: ", time.time() - t1)
    if params["do_ref"]:
        if step in params["reference_steps"]:
            print("Reference calculations for step", step)
            t1 = time.time()
            make_input(
                params["prefix"] + "_ref",
                params["reference_input_template"],
                params["reference_software"],
                params["trajectory_xyz_file"],
                step)
            # And run the calculations for this step again
            # if params["software"].lower()=="cp2k":
            run_single_point(params, params["prefix"] + "_ref", step, False)
            # os.system(F"mpirun -np {params['nprocs']} {params['software_exe']} -i input_{params['prefix']+'_ref'}_{step}.inp -o output_{params['prefix']+'_ref'}_{step}.log")
            save_data(params, params["prefix"] + '_ref', step, params["reference_dir"], False)
            print(f"Elapsed time for reference calculations: ", time.time() - t1)
    print(f"Done with step {step}!")


def find_steps(files):
    """
    This is an auxiliary function that finds how many steps were done
    in a directory by extracting the indices of the *.npy files
    Args:
        files (list): A list containing files in a give path - usually generated from glob library
    Returns:
        steps (list): The steps for which the calculations are complete and the data are present in a directory
    """
    steps = []
    for file in files:
        step = int(file.split('_')[-1].replace('.npy', ''))
        steps.append(step)
    steps = np.unique(steps)
    return steps


def distribute_jobs(params):
    """
    This function runs single-point electronic stucture calculations for
    geometries over a precomputed trajectory and distributes them over
    different ndoes
    Args:
        params (dictionary):
            prefix: The prefix used to names of the files when saved.
            trajectory_xyz_file: The full path to the precomputed trajectory `xyz` file.
            scratch: This flag is used to compute all the data, including both guess and reference data, from scratch.
            do_more: This flag is used to perform additional reference calculations. This may be needed when more data are required to train the ML model. This and `scratch` flags cannot have the same logical values at the same time.
            do_more_steps: The number of additional reference steps while the `do_more` flag is set to `True`.
            guess_dir: The full path to save the guess calculations.
            do_guess: A boolean flag for doing the guess calculations.
            do_ref: A boolean flag for doing only the reference calculations.
            reference_dir: The full path to save the reference calculations.
            guess_input_template: The input template used for perfroming guess calculations.
            reference_input_template: The input template required for performing reference calculations.
            guess_software: The software used to perform guess calculations. Current values are `cp2k` and `dftb+`
            guess_software_exe: The executable or the full executable path to the software required to compute the guess calculations.
            guess_mpi_exe: The `mpi` executable for running the `guess_software`.
            reference_software: The software to perform reference calculations. Current values are `cp2k` and `dftb+`.
            reference_software_exe: The executable or the full executable path to the software required to compute the reference calculations.
            reference_mpi_exe: The `mpi` executable for running the `reference_software`.
            reference_steps: The number of reference steps.
            njobs: The number of jobs for distributing the calculations.
            istep: The initial step in the precomputed molecular dynamics trajectory.
            fstep: The final step in the precomputed molecular dynamics trajectory.
            random_selection: A boolean flag for random shuffling of the istep to fstep.
            nprocs: The number of processors to be used for calculations in each job.
            remove_raw_outputs: This falg removes large raw files which are read by Libra and saved.
            submit_template: The full path to the submission file.
            submit_exe: The submission executable - for some HPC environments, it is `sbatch` and for some others, it is `qsub`. Currently, only slurm environment is tested.
            software_load_instructions: The loading instructions that are needed to load Python, guess and reference software, and any other thing required to run calculations.
    Returns:
        None
    """
    critical_params = ['guess_input_template', 'trajectory_xyz_file', 'reference_input_template', 'user_steps',
                       'submit_template', 'software_load_instructions']
    default_params = {
        'prefix': 'libra',
        'scratch': True,
        'do_more': False,
        'do_more_steps': 10,
        'random_selection': True,
        'guess_dir': './guess',
        'do_guess': True,
        'reference_dir': './ref',
        'do_ref': True,
        'guess_software': 'cp2k',
        'guess_software_exe': 'cp2k.psmp',
        'guess_mpi_exe': 'mpirun',
        'reference_software': 'cp2k',
        'reference_software_exe': 'cp2k.psmp',
        'reference_mpi_exe': 'mpirun',
        'reference_steps': 10,
        'njobs': 2,
        'nprocs': 2,
        'remove_raw_outputs': True,
        'submit_exe': 'sbatch'}
    # First load the default parameters etc
    comn.check_input(params, default_params, critical_params)
    # Now let's create the random numbers
    # nsteps = params["fstep"]-params["istep"]
    nsteps = len(params["user_steps"])
    params["istep"] = min(params["user_steps"])
    params["fstep"] = max(params["user_steps"])
    ref_steps = params["reference_steps"]
    if params["scratch"] == params["do_more"]:
        raise ("'scratch' and 'do_more' cannot get the same value at the same time!")
    if params["do_guess"] == False and params["do_ref"] == False:
        raise ("do_guess and do_ref cannot get False value at the same time!")
    # if params["scratch"]:
        # Remove evrything including the data
        # os.system("rm -rf job* sample_files")# shuffled_indices.npy") # Removes the job folders and the random indices
        # os.system('rm find . -name "*.npy" ') # This seems to be brutal! It removes everything :))
        # shuffled_indices = np.arange(params["istep"], params["fstep"])
        # params["steps"] = list(range(params["istep"], params["fstep"]))
        # if params["random_selection"]:
        #    np.random.shuffle(shuffled_indices)
        # np.save('shuffled_indices.npy', shuffled_indices)
        # params["reference_steps"] = shuffled_indices[0:ref_steps].tolist()
        # Here we add the first and last geometries
        # since this is an interpolation
        # if params["istep"] not in params["reference_steps"]:
        #    params["reference_steps"].append(params["istep"])
        # if params["fstep"] not in params["reference_steps"]:
        #    params["reference_steps"].append(params["fstep"])

    # elif params["do_more"]:
    #    # find how many steps were done in the reference_dir
    #    ref_files = glob.glob(f'{params["reference_dir"]}/*.npy')
    #    steps_done = find_steps(ref_files)
    #    #print(steps_done)
    #    if len(steps_done)==0:
    #        raise("Doing more steps is requested but the reference directory is empty!")
    #    ref_steps = len(steps_done)
    #    #print(ref_steps)
    #    more_steps = params["do_more_steps"]
    #    shuffled_indices = np.load("shuffled_indices.npy")
    #    params["reference_steps"] = shuffled_indices[ref_steps:ref_steps+more_steps].tolist()
    #    #print(params['reference_steps'])
    #    # For do_more, we don't need to add the first and last step since they are already computed in "scratch"
    # params["steps"] = np.sort(params["reference_steps"]).tolist() #
    # shuffled_indices[ref_steps:ref_steps+more_steps].tolist()

    # Create guess and reference directories
    os.system(f'mkdir {params["reference_dir"]} {params["guess_dir"]}')
    if len(params["user_steps"]) > 0:
        params["steps"] = params["user_steps"]
        params["reference_steps"] = params["user_steps"]
    if (len(params["steps"]) / params["njobs"]) < 2:
        raise ("The division of steps/njobs is less than 2! Please select smaller number of jobs.")
    steps_split = np.array_split(params["steps"], params["njobs"])
    # os.system(F"mkdir {res_directory}")
    # Read the submit file
    f = open(params["submit_template"], "r")
    lines_submit = f.readlines()
    f.close()

    for job in range(params["njobs"]):
        print(F"Submitting job {job}...")
        if params["scratch"]:
            os.system(F"mkdir job_{job}")
            os.chdir(F"job_{job}")
        elif params["do_more"]:
            os.system(F"mkdir job_more_{job}")
            os.chdir(F"job_more_{job}")
        this_job_steps = steps_split[job].tolist()
        # params["job_istep"] = int(this_job_steps[0])
        # params["job_fstep"] = int(this_job_steps[-1])
        params["job_steps"] = this_job_steps
        # Now dump the params information
        with open("params.json", "w") as f:
            #    print(params)
            json.dump(params, f)
        # np.savetxt("steps.txt", this_job_steps, fmt="%d")
        f = open(F"submit.slm", "w")
        for i in range(len(lines_submit)):
            if "#" in lines_submit[i]:
                f.write(lines_submit[i])
        f.write("\n\n\n")
        f.write(params["software_load_instructions"])
        f.write("\n")
        f.write("python run.py \n\n")
        f.close()
        # Now the python file: run.py
        f = open(f"run.py", "w")
        f.write("""import json
from libra_py.workflows.nbra.generate_data import *
with open('params.json', 'r') as f:
    params = json.load(f)
for step in params['job_steps']:
    gen_data(params, step)
        """)
        f.close()
        os.system(F"{params['submit_exe']} submit.slm")
        os.chdir("../")
