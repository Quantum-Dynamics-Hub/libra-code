#***********************************************************
# * Copyright (C) 2013-2019 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/
"""
.. module:: hpc_utils
   :platform: Unix, Windows
   :synopsis: 
       This module implements the functionality to run distributed calculations
       on the HPC

.. moduleauthor:: Alexey V. Akimov

"""


import os
import sys
import time
import re


def substitute(filename_template, out_filename, params):
    """

    Make a copy of a template file but with some variables replaced

    Args:
        filename_template ( string ): the name of the file that serves as a template
        out_filename ( string ): the name of the output (result) file
        params ( dictionary ): the keys of this dictionary should are the strings 
            (more generally regular expressions) looked for, the values - are the strings
            with which found regexes will be substituted.

    Returns:
        None: but creates the out_filename files

    """

    f = open(filename_template,"r")
    A = f.readlines()
    f.close()

    f_out = open("%s" % (out_filename), "w")

    for a in A:        
        b = a
        for key in params.keys():
            b = re.sub(key, params[key], b)

        f_out.write(b)
    f_out.close()



def make_submit(Nstart,Nend,job_dir,submit_templ):
    """    
  
    This function makes a PBS/SLURM submit file in a given job directory
    using a template file. The function looks for the placeholders 
    specified by "param1= " and "param2= " and appends actual values
    of these parameters as needed.

    Args: 
        Nstart ( int ): index of a starting parameter 
        Nend ( int ): index of a finishing parameter 
        job_dir ( string ): name of the directory where this file will be printed out to       
        submit_templ ( string ): the template file name

    Returns:
        None: but generates actual submit files

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

    This function prepares a group of input files to be executed as a single job
    note the <job> here refers to a single PBS/SLURM submission file and may be executed
    on several processors.

    Args: 
        Nstart ( int ): index of a starting parameter 
        Nend ( int ): index of a finishing parameter 
        job_dir ( string ): name of the directory where this file will be printed out to       
        prefix ( string ): the template file name [e.g. "x.scf"].

    Returns:
        None: just moves the files

    Example:    
        >>> job(0, 5, "job0", "x.scf")
        Will:
        1. Create a directory "job0"
        2. Move filex x.scf.0.in to x.scf.5.in to the directory "job0"

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

    This function creates a number of jobs from the pool of input files 
    which start from ```Nmin``` to ```Nmax``` with the maximal number 
    of input files in one job given by ```max_steps```.  
    It also starts the job in a given directory

    Args: 
        Nmin ( int ): minimal index of the file (a part of the file name) in the pool
        Nmax ( int ): maximal index of the file (a part of the file name) in the pool
        max_steps ( int ): how many files to handle in one job
        submit_templ ( string ): the filename of the PBS/SLURM submit file that 
            contains all the information on how to run the jobs involving a subset of
            files, but doesn't specify which particular files to handle. These parameters
            will be automatically setup by this funciton           
        exp_files ( list of strings ): is a list of the input files defining how to do an 
            export of QE wavefunctions
        prefixes ( list of strings ): is a list of prefixes of the files to be distributed
        do_submit ( int ): a flag to choose if we actually want to submit the jobs
            (do_submit==1 || ==2 || ==3 ) or only distribute the files (otherwise):

            - 0: only discribute files, no actual execution
            - 1: distribute and submit using PBS
            - 2: distribute and submit using SLURM
            - 3: distribute and submit using PYTHON no scheduling systems - good for runs on a local computer

    Returns:
        None: just organizes the execution of the calculations 

    """

    j = 0  # job index
    Nstart = 0
    Nend = max_steps-1
    njobs = int( (Nmax -1 - Nmin)/max_steps )

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


