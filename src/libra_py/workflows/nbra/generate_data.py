import os
import numpy as np
import sys
import libra_py.packages.cp2k.methods as CP2K_methods

def make_input(prefix, input_template, software, trajectory_xyz_file, step):

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


def run_calculations(params):

    prefix = params["prefix"] 
    trajectory_xyz_file = params["trajectory_xyz_file"]
    input_template = params["input_template"] 
    software = params["software"] 
    steps = params["steps"] 
    software_exe = params["software_exe"] 
    nprocs = params["nprocs"]
    remove_raw_outputs = params["remove_raw_outputs"]
    res_directory = params["res_directory"]
    for step in steps:
        make_input(prefix, input_template, software, trajectory_xyz_file, step)
        if software.lower()=="cp2k":
            os.system("export OMP_NUM_THREADS=1")
            os.system(F"mpirun -np {nprocs} {software_exe} -i input_{prefix}_{step}.inp -o output_{prefix}_{step}.log")
            try:
                filename = F"{prefix}_{step}-libra-1_0.Log"
                properties, data = CP2K_methods.read_ao_matrices(filename)
                for i in range(len(data)):
                    tmp_name = "".join(properties[i].split()).replace("-","_").lower()
                    np.save(F"{res_directory}/{prefix}_{tmp_name}_{step}.npy", data[i])
                if remove_raw_outputs:
                    os.system(F"rm {prefix}_{step}-libra-*.Log {prefix}_{step}-RESTART.wfn*")
            except:
                raise("Could not read the data from CP2K .Log files!")
        elif software.lower()=="dftb+":
            os.system(F"export OMP_NUM_THREADS={nprocs}") 
            os.system(F"mv dftb_in_{prefix}_{step}.hsd dftb_in.hsd") 
            os.system(F"{software_exe} dftb_in.hsd > output_{prefix}_{i}.log") 
            os.system(F"mv dftb_in.hsd dftb_in_{prefix}_{step}.hsd") 
            try:
                filename = "hamsqrd1.dat"
                if os.path.exists(filename):
                    data = np.loadtxt(filename, skiprows=5)
                    np.save(F"{prefix}_ham_{step}.npy", data)
                else:
                    raise("Could not find file hamsqrd1.dat...")
                #filename = "hamsqrd1.dat"
                #if os.path.exists(filename):
                #    data = np.loadtxt(filename, skiprows=5)
                #    np.save("{prefix}_ham_{step}.npy", data)
            except:
                raise("Could not read data from DFTB+ .dat files!")

        
def distribute_jobs(params):
    steps = params["steps"] 
    njobs = params["njobs"]
    submit_template = params["submit_template"]
    run_template = params["run_template"]
    software_load_instructions = params["software_load_instructions"]
    submit_exe = params["submit_exe"] 
    res_directory = params["res_directory"]
    f = open(submit_template, "r")
    lines_submit = f.readlines()
    f.close()

    steps_split = np.array_split(steps, njobs)
    os.system(F"mkdir {res_directory}")
    for job in range(njobs):
        print(F"Submitting job {job}...")
        os.system(F"mkdir job_{job}")
        os.chdir(F"job_{job}")
        this_job_steps = steps_split[job]
        np.savetxt("steps.txt", this_job_steps, fmt="%d")
        f = open(F"submit.slm","w")
        for i in range(len(lines_submit)):
            if "#" in lines_submit[i]:
                f.write(lines_submit[i])
        f.write("\n\n\n")
        f.write(software_load_instructions)
        f.write("\n")
        f.write("python run.py \n\n")
        os.system(F"cp ../{run_template} run.py")
        f.close()
        os.system(F"{submit_exe} submit.slm")
        os.chdir("../")



