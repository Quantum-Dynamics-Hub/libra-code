#***********************************************************
# * Copyright (C) 2024 Mohammad Shakiba and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

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
    """
    f = open(input_template,'r')
    lines_input = f.readlines()
    f.close()

    CP2K_methods.read_trajectory_xyz_file(trajectory_xyz_file, step)

    if software.lower()=="cp2k":
        f = open(F"input_{prefix}_{step}.inp", "w")
        for i in range(len(lines_input)):
            if "COORD_FILE_NAME".lower() in lines_input[i].lower():
                f.write(F"     COORD_FILE_NAME  coord-{step}.xyz \n")
            elif "project" in lines_input[i].lower() or "project_name" in lines_input[i].lower():
                f.write(F"  PROJECT {prefix}_{step}\n")
            elif "filename" in lines_input[i].lower() and "!" not in lines_input[i]:
                f.write("     FILENAME libra\n")
            else:
                f.write(lines_input[i])
        f.close()
    elif software.lower()=="dftb+":
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
        for i in range(end+1):
            f.write(lines_input[i])
        f.close()

def save_data(params, prefix, step, directory):
    if params["software"].lower()=="cp2k":
        try:
            filename = F"{prefix}_{step}-libra-1_0.Log"
            properties, data = CP2K_methods.read_ao_matrices(filename)
            for i in range(len(data)):
                tmp_name = "".join(properties[i].split()).replace("-","_").lower()
                np.save(F"{directory}/{prefix}_{tmp_name}_{step}.npy", data[i])
            if step==params["istep"]: 
                # This is how we can get the sample molden, wfn, and Log files.
                os.system("mkdir ../sample_files")
                os.system("mv {prefix}_{step}* ../sample_files")
            if params["remove_raw_outputs"]: # and step>params["istep"]:
                os.system(F"rm {prefix}_{step}-libra-*.Log {prefix}_{step}-RESTART.wfn* {prefix}_{step}*molden")
        except:
            raise("Could not read the data from CP2K .Log files!")

    elif params["software"].lower()=="dftb+":
        try:
            filename = "hamsqrd1.dat"
            if os.path.exists(filename):
                data = np.loadtxt(filename, skiprows=5)
                np.save(F"{prefix}_ham_{step}.npy", data)
            else:
                raise("Could not find file hamsqrd1.dat...")
        except:
            raise("Could not read data from DFTB+ .dat files!")

def run_single_point(params, prefix, step):
    """
    This function runs the software on the Linux environment
    """
    if params["software"].lower()=="cp2k":
        os.system(F"mpirun -np {params['nprocs']} {params['software_exe']} -i input_{prefix}_{step}.inp -o output_{prefix}_{step}.log")
    elif params["software"].lower()=="dftb+":
        os.system(F"export OMP_NUM_THREADS={nprocs}") 
        os.system(F"mv dftb_in_{prefix}_{step}.hsd dftb_in.hsd") 
        os.system(F"{software_exe} dftb_in.hsd > output_{prefix}_{i}.log") 
        os.system(F"mv dftb_in.hsd dftb_in_{prefix}_{step}.hsd") 


def gen_data(params, step):
    """
    This function takes multiple parameters to run electronic structure calculations
    """
    # Generate the step^th geometry from the trajectory file
    # and create an input file for that
    make_input(params["prefix"], params["guess_input_template"], params["software"], params["trajectory_xyz_file"], step)
    # Run the calculations for this input file
    t1 = time.time()
    print("Guess calculations for step", step)
    run_single_point(params, params["prefix"], step)

    save_data(params, params["prefix"], step, params["guess_dir"])
    #else:
    #    print("to be implemented!")
    print(f"Elapsed time for guess calculations: ", time.time()-t1)
    # After the calculation is done for the guess geometry
    # if the step is in the reference steps, create another input file
    if step in params["reference_steps"]:
        print("Reference calculations for step", step)
        t1 = time.time()
        make_input(params["prefix"]+"_ref", params["reference_input_template"], params["software"], params["trajectory_xyz_file"], step)
        # And run the calculations for this step again
        #if params["software"].lower()=="cp2k":
        run_single_point(params, params["prefix"]+"_ref", step)
        #os.system(F"mpirun -np {params['nprocs']} {params['software_exe']} -i input_{params['prefix']+'_ref'}_{step}.inp -o output_{params['prefix']+'_ref'}_{step}.log")
        save_data(params, params["prefix"]+'_ref', step, params["reference_dir"])
        print(f"Elapsed time for reference calculations: ", time.time()-t1)
    print(f"Done with step {step}!")


def find_steps(files):
    """
    This is an auxiliary function that finds how many steps were done
    in a directory by extracting the indices of the *.npy files
    """
    steps = []
    for file in files:
        step = int(file.split('_')[-1].replace('.npy',''))
        steps.append(step)
    steps = np.unique(steps)
    return steps
    
def distribute_jobs(params):
    """
    This function runs single-point electronic stucture calculations for 
    geometries over a precomputed trajectory and distributes them over 
    different ndoes
    """
#    critical_params = ['lowest_orbital','highest_orbital']
#    default_params = {
#                     }
    # First load the default parameters etc
#    comn.check_input(params, default_params, critical_params)
    # Now let's create the random numbers
    nsteps = params["fstep"]-params["istep"]
    ref_steps = params["reference_steps"]
    if params["scratch"]:
        # Remove evrything including the data
        os.system("rm -rf job*") # Removes the job folders
        #os.system('rm find . -name "*.npy" ') # This seems to be brutal! It removes everything :))
        shuffled_indices = np.arange(params["istep"], params["fstep"])
        params["steps"] = list(range(params["istep"], params["fstep"]))
        np.random.shuffle(shuffled_indices)
        np.save('shuffled_indices.npy', shuffled_indices)
        params["reference_steps"] = shuffled_indices[0:ref_steps].tolist()
        # Here we add the first and last geometries 
        # since this is an interpolation
        if params["istep"] not in params["reference_steps"]:
            # This is how it is done in numpy
            params["reference_steps"] = np.append(params["reference_steps"], params["istep"])
        if params["fstep"] not in params["reference_steps"]:
            params["reference_steps"] = np.append(params["reference_steps"], params["fstep"])

    elif params["do_more"]:
        # find how many steps were done in the reference_dir
        ref_files = glob.glob(f'{params["reference_dir"]}/*.npy')
        steps_done = find_steps(ref_files)
        if steps_done==0:
            raise("Doing more steps is requested but the reference directory is empty!")
        ref_steps = steps_done
        more_steps = params["do_more_steps"]
        shuffled_indices = np.load("shuffled_indices.npy")
        params["reference_steps"] = shuffled_indices[ref_steps:ref_steps+more_steps].tolist()
        # For do_more, we don't need to add the first and last step since they are already computed in "scratch"
        params["steps"] = shuffled_indices[ref_steps:ref_steps+more_steps].tolist()
        
    # Create guess and reference directories
    os.system(f'mkdir {params["reference_dir"]} {params["guess_dir"]}')
    steps_split = np.array_split(params["steps"], params["njobs"])
    #os.system(F"mkdir {res_directory}")
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
        params["job_istep"] = int(this_job_steps[0])
        params["job_fstep"] = int(this_job_steps[-1])
        # Now dump the params information
        with open("params.json", "w") as f:
        #    print(params)
            json.dump(params, f)
        #np.savetxt("steps.txt", this_job_steps, fmt="%d")
        f = open(F"submit.slm","w")
        for i in range(len(lines_submit)):
            if "#" in lines_submit[i]:
                f.write(lines_submit[i])
        f.write("\n\n\n")
        f.write(params["software_load_instructions"])
        f.write("\n")
        f.write("python run.py \n\n")
        f.close()
        # Now the python file: run.py
        f = open(f"run.py","w")
        f.write("""import json
from libra_py.workflows.nbra.generate_data import *
with open('params.json', 'r') as f:
    params = json.load(f)
for step in range(params['job_istep'], params['job_fstep']):
    gen_data(params, step)
        """)
        f.close()
        os.system(F"{submit_exe} submit.slm")
        os.chdir("../")

