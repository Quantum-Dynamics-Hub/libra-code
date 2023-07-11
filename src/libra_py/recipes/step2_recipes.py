#***********************************************************
# * Copyright (C) 2021 Brendan Smith and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys
import glob
import logging
import numpy as np
import multiprocessing as mp

import matplotlib.pyplot as plt

import util.libutil as comn
import libra_py.packages.cp2k.methods as CP2K_methods
import libra_py.packages.gaussian.methods as Gaussian_methods
import libra_py.packages.dftbplus.methods as DFTB_methods
from libra_py.workflows.nbra import step2_many_body


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler('step2_recipes.log')
file_handler.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def make_step2_submit_template(params):
    """ Generate the submit_template.slm file for step2 in the nbra workflow.
            
    Please note that the generated submit_template.slm here is
    to be used on the UB-CCR.   

    Args:
        params (dict): The parameters dictionary used for containing 
            the variables needed to make the submit.slm template 

    Returns:
        None: but creates a file template for submit.slm.

    Notes:
        cp2k_exe=/panasas/scratch/grp-alexeyak/brendan/cp2k/exe/local/cp2k.popt
        gaussian_exe=g16
        dftb_exe=/util/academic/dftbplus/19.1-arpack/bin/dftb+
 
    """

    logger.debug("Entered into generate_step2_submit_template function") 

    critical_params = []
    default_params = { "partition":"debug", "clusters":"ub-hpc", "constraint":"CPU-Gold-6130",
                       "nnodes":1, "ncores":1, "cpus_per_task":1,
                       "time":"1:00:00", "memory":50000, "email":None, 

                       "es_software":"cp2k",  "es_software_input_template": "cp2k_input_template.inp",
                       "es_software_exe":"/panasas/scratch/grp-alexeyak/brendan/cp2k/exe/local/cp2k.popt",
                       "waveplot_exe":"/util/academic/dftbplus/20.2.1-arpack/bin/waveplot",

                       "project_name":"libra_step2",
                       "min_band":0, "max_band":1, "homo_index":0, 
                       "isUKS":0,   

                       "trajectory_xyz_filename":"md.xyz", "results_path":"",

                       "do_cube_visualization":0, "states_to_be_plotted":0
                      }
    comn.check_input(params, default_params, critical_params)
    logger.debug("Checked params in function generate_template")

    file_content = "#!/bin/sh\n#SBATCH --requeue\n"

    #============ General Slurm stuff ===============
    if params["partition"]!=None:
        partition = params["partition"]
        file_content = file_content + F"""#SBATCH --partition={partition} --qos={partition}\n"""

    if params["clusters"]!=None:
        clusters = params["clusters"]
        file_content = file_content + F"""#SBATCH --clusters={clusters}\n"""

    if params["constraint"]!=None:
        constraint = params["constraint"]
        file_content = file_content + F"""#SBATCH --constraint={constraint}\n"""

    if params["nnodes"]!=None:
        nnodes = params["nnodes"]
        file_content = file_content + F"""#SBATCH --nodes={nnodes}\n"""

    if params["ncores"]!=None:
        ncores = params["ncores"]
        file_content = file_content + F"""#SBATCH --ntasks-per-node={ncores}\n"""

    if params["cpus_per_task"]!=None:
        cpus_per_task = params["cpus_per_task"]
        file_content = file_content + F"""#SBATCH --cpus-per-task={cpus_per_task}\n"""

    if params["time"]!=None:
        time = params["time"]
        file_content = file_content + F"""#SBATCH --time={time}\n"""

    if params["memory"]!=None:
        memory = params["memory"]
        file_content = file_content + F"""#SBATCH --mem={memory}\n"""

    if params["email"]!=None:
        email = params["email"]
        file_content = file_content + F"""#SBATCH --mail-user={email}\n"""


    #========== Additional setup for CCR  ==================

    file_content = file_content + F"""
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST="$SLURM_JOB_NODELIST
echo "SLURM_NNODES="$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory="$SLURM_SUBMIT_DIR

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads

"""

    #===== Module loading ===========
    modules = params["modules"]
    for module in modules:
        file_content = file_content + F"""module load {module}\n"""

    #===== Initialize variables to be updated by the job distributor ===========
    file_content = file_content + F"""

job_init_step=
nsteps_this_job=
njob=

"""

    #================= Python script to run ==========================
    es_software = params["es_software"]
    es_software_input_template = params["es_software_input_template"]
    es_software_exe = params["es_software_exe"]
    project_name = params["project_name"]
    trajectory_xyz_filename = params["trajectory_xyz_filename"]

    isUKS = params["isUKS"]
    min_band = params["min_band"]
    max_band = params["max_band"]
    homo_index = params["homo_index"]

    path = params["results_path"]

    do_cube_visualization = params["do_cube_visualization"]
    states_to_be_plotted = params["states_to_be_plotted"]

    file_content = file_content +F"""
python -c "from libra_py.workflows.nbra import step2_many_body
params = {{}}
params[\\"es_software\\"]=\\"{es_software}\\"
params[\\"es_software_input_template\\"]=\\"{es_software_input_template}\\"
params[\\"es_software_exe\\"]=\\"{es_software_exe}\\" """

    if es_software == "dftb+":    
        waveplot_exe = params["waveplot_exe"]
        file_content = file_content +F"""
params[\\"waveplot_exe\\"]=\\"{waveplot_exe}\\" """

    file_content = file_content +F"""
params[\\"nprocs\\"]=\\"$SLURM_NTASKS_PER_NODE\\"
params[\\"project_name\\"]=\\"{project_name}\\"
params[\\"trajectory_xyz_filename\\"]=\\"{trajectory_xyz_filename}\\"
params[\\"isUKS\\"]={isUKS}
params[\\"min_band\\"]={min_band}
params[\\"max_band\\"]={max_band}
params[\\"ks_orbital_homo_index\\"]={homo_index}
params[\\"istep\\"]=\\"$job_init_step\\"
params[\\"nsteps_this_job\\"]=\\"$nsteps_this_job\\"
params[\\"njob\\"]=\\"$njob\\"
params[\\"res_dir\\"]=\\"{path}/res\\"
params[\\"do_cube_visualization\\"]={do_cube_visualization}
params[\\"path_to_tcl_file\\"] = \\"{path}/cube.tcl\\"
params[\\"states_to_be_plotted\\"]=\\"{states_to_be_plotted}\\"
params[\\"mo_images_directory\\"] = \\"{path}/mo_images\\"
params[\\"logfile_directory\\"]=\\"logfiles\\"
print( params )
step2_many_body.run_step2_many_body( params )
"

"""   

    f = open("submit_template.slm", "w")
    f.write(file_content)
    f.close()


def run_step2_jobs(params):
    """ Prepares and runs job folders needed for running step2 

    Args:
        params (dict): The parameters dictionary used for containing 
            the variables needed to make the submit.slm template 

    Returns:
        None: but creates the wd (working diretory) that stores the 
            job folders and the files needed to run each job, and runs
            those jobs.    
    """

    logger.debug("Entered into the function initialize_step2_jobs")

    critical_params = []
    default_params = { "trajectory_xyz_filename":"md.xyz", "es_software":"cp2k", "es_software_input_template":"cp2k_input_template", "istep":0, "fstep":1, "njobs":1, "waveplot_input_template":"waveplot_in.hsd" }
    comn.check_input(params, default_params, critical_params)
    logger.debug("Checked params in the function initialize_step2_jobs")
    
    # Consider using shutils module here instead of os.system()
    os.system("rm -r wd")
    os.mkdir('wd')
    os.system("rm -r res")
    os.mkdir('res')

    # Get the neededvariables from the params dictionary
    trajectory_xyz_filename = params["trajectory_xyz_filename"]
    es_software = params["es_software"]
    es_software_input_template = params["es_software_input_template"]
    init_md_step = params["istep"]
    final_md_step = params["fstep"]
    njobs = params["njobs"]

    if es_software == "dftb+":
        waveplot_input_template = params["waveplot_input_template"]

    # Initialize the jobs
    for njob in range(njobs):

        logger.debug(f"Within the loops over njobs, initializing job: {njob}")

        job_init_step, job_final_step = step2_many_body.curr_and_final_step_job( init_md_step, final_md_step, njobs, njob )
        nsteps_this_job = job_final_step - job_init_step + 1
        logger.debug(f"The initial step for the job {njob} is: {job_init_step} with final step: {job_final_step}")
        logger.debug(f"nsteps_this_job is: {nsteps_this_job}")

        logger.debug(f"Creating job: {njob} in directory: {os.getcwd()}")
        if es_software.lower() == "cp2k":
            logger.debug("es_software.lower() == cp2k")
            CP2K_methods.cp2k_distribute( job_init_step, job_final_step, nsteps_this_job, trajectory_xyz_filename, es_software_input_template, njob )

        elif es_software.lower() == "gaussian":
            logger.debug("es_software.lower() == gaussian")
            Gaussian_methods.gaussian_distribute( job_init_step, job_final_step, nsteps_this_job, trajectory_xyz_filename, es_software_input_template, njob )

        elif es_software.lower() == "dftb+":
            logger.debug("es_software.lower() == dftb+")
            DFTB_methods.dftb_distribute( job_init_step, job_final_step, nsteps_this_job, trajectory_xyz_filename, es_software_input_template, waveplot_input_template, njob )

        logger.debug(f"Finished distributing for job: {njob}. We are currently in directory: {os.getcwd()}")
        os.chdir("wd/job"+str(njob)+"/")
        logger.debug(f"Changed directory to: {os.getcwd()}")


        logger.debug("Assigning values to the variables job_init_step, nsteps_this_job, and njob in submit_template.slm")

        os.system("cp ../../submit_template.slm submit_"+str(njob)+".slm")
        # Now, open the submit_template file in this job folder
        # Add the values to the params for this job
        f = open( "submit_"+str(njob)+".slm" )
        submit_template_file = f.readlines()
        submit_template_file_size = len( submit_template_file )
        f.close()

        f = open( "submit_"+str(njob)+".slm" , 'w' ); f.close()
        for i in range( submit_template_file_size ):      

            f = open( "submit_"+str(njob)+".slm" , 'a' )

            submit_template_file_line = submit_template_file[i].split()
            if not submit_template_file_line:
                continue

            elif submit_template_file_line[0] == "job_init_step=":
                f.write( "declare -i job_init_step=%i" % (job_init_step) )
                f.write("\n")

            elif submit_template_file_line[0] == "nsteps_this_job=":
                f.write( "declare -i nsteps_this_job=%i" % (nsteps_this_job) )
                f.write("\n")

            elif submit_template_file_line[0] == "njob=":
                f.write( "declare -i njob=%i" % (njob) )
                f.write("\n")

            else:
                f.write( submit_template_file[i] )            

        logger.debug(f"Updated submit_template.slm")
 
        os.system("sbatch submit_"+str(njob)+".slm")
        #os.system("sh submit_"+str(njob)+".slm")
        logger.debug(f"Submitting the job in folder: {os.getcwd()}")

        # Change directory back to the top directory
        os.chdir("../../")
        logger.debug(f"Finished submitting job. We are now back in directory: {os.getcwd()}")



if __name__ == '__main__':
    params = { "partition":"debug", "clusters":"ub-hpc", "time":"1:00:00", "ncores":32, "memory":50000, "mail":"bsmith24@buffalo.edu", "min_band":28, "max_band":29, "homo_index":28, "project_name":"c10h16", "trajectory_xyz_filename":"c10h16.xyz", "path":"/panasas/scratch/grp-alexeyak/brendan/active_projects/libra_development/test_black_box", "es_software":"cp2k", "es_software_input_template":"cp2k_input_template.inp", "istep":0, "fstep":1, "njobs":1}
    logger.debug("Running step2 recipe")

    make_step2_submit_template(params)
    run_step2_jobs(params)
