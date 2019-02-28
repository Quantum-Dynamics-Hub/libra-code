#***********************************************************
# * Copyright (C) 2013-2018 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys
import time

def make_submit(Nstart,Nend,job_dir,submit_templ):
    """
    This function makes the submit file in the job_dir
    submit_templ - is a file name of the template for submit file
    """

    f = open(submit_templ,"r")
    A = f.readlines()
    f.close()

    f_out = open("%s/%s" % (job_dir,submit_templ), "w")
    for a in A:
        if a.find("param1=")!=-1:
            a = a[:-1] + str(Nstart) + "\n"
        elif a.find("param2=")!=-1:
            a = a[:-1] + str(Nend) + "\n"
        f_out.write(a)
    f_out.write("\n")
    f_out.close();


def job(Nstart,Nend,job_dir,prefixes):
    """
    This function prepares the group of input files to be executed as a single job
    note the <job> here referres to a single PBS submission file and may be executed
    on several processors
    Nstart - is the minimal index of the job file
    Nend - is the maximal index of the job file
    prefix - is input file prefix (was x.scf)
    """
    # Create <job_dir> directory if it does not exist yet
    if os.path.isdir(job_dir):
        pass
    else:
        os.system("mkdir %s" % job_dir)

    # Copy all files
    i = Nstart
    while i<=Nend:
        for prefix in prefixes:
            os.system("cp %s.%d.in %s" % (prefix,i,job_dir))
        i = i + 1


def distribute(Nmin,Nmax,max_steps,submit_templ,exp_files,prefixes,do_submit):
    """
    This function creates a number of jobs from the pool
    of the input files which start from Nmin to Nmax with the
    maximal number of input files in one job given by max_steps
    It also starts the job in a given directory
    exp_files - is a list of the input files for export
    prefixes - is a list of prefixes of the files to be distributed
    do_submit - a flag to choose if we actually want to submit the jobs (do_submit==1 || ==2)
                or only distribute the files (otherwise)
    do_submit == 0 - only discribute files, no actual execution
    do_submit == 1 - distribute and submit using PBS
    do_submit == 2 - distribute and submit using SLURM
    do_submit == 3 - distribute and submit using PYTHON no scheduling systems - good for runs on a local computer

    """

    j = 0  # job index
    Nstart = 0
    Nend = max_steps-1
    njobs = (Nmax -1 - Nmin)/max_steps

    while j<njobs:
        job(Nstart,Nend+1,"job%d" % j,prefixes) # add 1 to avoid merging adjacent runs
        if do_submit in [0, 1, 2]:
            make_submit(Nstart,Nend+1,"job%d" % j,submit_templ)

        Nstart = Nstart + max_steps 
        Nend = Nend + max_steps    
        # Go into that directory and submit the job  
        os.chdir("job%d" % j)
       
        for exp_file in exp_files:
            os.system("cp ../%s ." % exp_file)

        if do_submit==1:
            os.system("qsub %s" % submit_templ)
            time.sleep(10)
        elif do_submit==2:
            os.system("sbatch %s" % submit_templ)
            time.sleep(10)

        os.chdir("../")
        j = j + 1

    job(Nstart,min(Nend,Nmax),"job%d" % j,prefixes)

    if do_submit in [0, 1, 2]:
        make_submit(Nstart,min(Nend,Nmax),"job%d" % j,submit_templ)
    os.chdir("job%d" % j)

    for exp_file in exp_files:
        os.system("cp ../%s ." % exp_file)

    if do_submit==1:
        os.system("qsub %s" % submit_templ)    
        time.sleep(10)  # we need to wait some time before submitting a new job - to make sure the memory is available
    elif do_submit==2:
        os.system("sbatch %s" % submit_templ)
        time.sleep(10)

    os.chdir("../")


