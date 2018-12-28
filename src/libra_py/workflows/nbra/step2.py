#***********************************************************
# * Copyright (C) 2017-2018 Brendan A. Smith, Wei Li, and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

from utils import *
import common_utils as comn
import compute_properties
import compute_hprime


def run_qe(params, t, dirname0, dirname1):
    """
    This fucntions runs the QE code for a system whose information is stored in the 
    params dictionary

    \param[in] params A dictionary containing important simulation parameters
    \param[in] t The current time step
    \param[in] dirname0 Name of the temporary directory where data will be stored 
              for the case without SOC 
    \param[in] dirname1 Name of the temporary directory where data will be stored 
              for the case with SOC 
    """

    tim = Timer()
    tim.start()

    # Now try to get parameters from the input
    BATCH_SYSTEM = get_value(params,"BATCH_SYSTEM","srun","s")  # either "srun" (for SLURM) or "mpirun" (for PBS)
    NP = get_value(params,"NP","1","i")
    EXE = get_value(params,"EXE","","s")
    EXE_EXPORT = get_value(params,"EXE_EXPORT","","s")
    prefix0 = get_value(params,"prefix0","x0.scf","s")
    prefix1 = get_value(params,"prefix1","x1.scf","s")
    nac_method = get_value(params,"nac_method",0,"i")  # choose what method for NAC calculations to use: 0 -standard, 1-corrected
    wd = get_value(params,"wd","wd","s")

    # Run calculations
    # A regular calculation anyway
    if nac_method == 0 or nac_method == 1 or nac_method == 3:
        os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix0,t,prefix0,t) )
        os.system( "%s -n %s %s < x0.exp.in > x0.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

        # Create temporary directory
        os.system("mkdir %s/%s" % (wd, dirname0) )

        # Copy some results to that directory
        os.system( "mv %s.%d.out %s/%s" % (prefix0,t, wd, dirname0) )
        os.system( "mv *.wfc* %s/%s" % (wd, dirname0) )
        os.system( "mv x0.export %s/%s" % (wd, dirname0) ) # "x0" - corresponds to x0 as a prefix in input files
                                                                        
    # Perform the soc calculation on its own, or in addition to the regular one
    if nac_method == 2 or nac_method == 3:
        os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix1,t,prefix1,t) )
        os.system( "%s -n %s %s < x1.exp.in > x1.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

        os.system("mkdir %s/%s" % (wd,dirname1) )

        os.system( "mv %s.%d.out %s/%s" % (prefix1,t, wd, dirname1) )
        os.system( "mv *.wfc* %s/%s" % (wd, dirname1) )
        os.system( "mv x1.export %s/%s" % (wd, dirname1) ) # "x1" - corresponds to x1 as a prefix in input files

    print "The time to run the QE calculations = ", tim.stop(); 



def read_info(params):
    """
    This fucntions reads the output from a QE calculations, and stores the output
    information in dictionaries

    \param[in] params A dictionary containing important simulation parameters
  
    Returns: Lists that contains system specific information such as the values
             for lattice constants, eigenvalues, etc. 
    """

    tim = Timer()
    tim.start()

    nac_method = get_value(params,"nac_method",0,"i")  # choose what method for NAC calculations to use: 0 -standard, 1-corrected
    wd0 = get_value(params,"wd0","wd","s")
    wd1 = get_value(params,"wd1","wd","s")

    info0, all_e_dum0 = None, None
    info1, all_e_dum1 = None, None

    # for non-relativistic, non-spin-polarized calculations
    if nac_method == 0:
        info0, all_e_dum0 = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd0, [], 0)
#        info0, all_e_dum0 = QE_methods.read_qe_index("%s/index.xml" % wd0, [], 0)

        if info0["nspin"] != 1:
            print "Error, you are not running the non-relativistic, non-spin-polarized calculation \
                   check your setting with nspin"
            sys.exit(0)
        print "The total # of k-points (non-spin-polarized calculation) is: ", info0["nk"]

    # for non-relativistic, spin-polarized calculations
    if nac_method == 1 or nac_method == 3:
        info0, all_e_dum0 = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd0, [], 0)
#        info0, all_e_dum0 = QE_methods.read_qe_index("%s/index.xml" % wd0, [], 0)

        if info0["nspin"] != 2:
            print "Error, you are not running spin-polarized calculations,\
                   check you settings with nspin"
            sys.exit(0)
        print "The total # of k-points (spin-polarized) including up and down components is: ", info0["nk"]

    # for fully-relativistic non-collinear calculations
    if nac_method == 2 or nac_method==3:
        info1, all_e_dum1 = QE_methods.read_qe_index("%s/curr1/x1.export/index.xml" % wd0, [], 0)
#        info1, all_e_dum1 = QE_methods.read_qe_index("%s/index.xml" % wd1, [], 0)

        if info1["nspin"] != 4:
            print "Error,you are not running SOC calculations \
                   check you setting with nspin. Also, veriy that you use fully-relativistic pseudopotentials"
            sys.exit(0)
        print "The total # of k-points (soc) is: ", info1["nk"]

    print "The time to get the basic parameters about your QE calculations = ", tim.stop();  

    return info0, all_e_dum0, info1, all_e_dum1
    



def read_all(params):
    """
    This function reads index, wfc and grid files from a given directory
    The number of wfc and grid files may be larger than 1 - this is the
    case of spin-polarized or multiple k-points calculations

   \param[in] params A dictionary containing important simulation parameters
              params["prefix"] - (string) - the location of the folder containing index.xml, wfc.*, and grid.* files
              params["read_wfc"] - (0 or 1) - whether or not to read the wfc coefficients. Default: 1
              params["read_grid"] - (0 or 1) - whether or not to read the grid informations. Default: 1
              params["verb0"] - (0 or 1) - turn off/on the extra printout while reading index.xml. Default: 0
              params["verb1"] - (0 or 1) - turn off/on the extra printout while reading wfc.*. Default: 0
              params["verb2"] - (0 or 1) - turn off/on the extra printout while reading grid.*. Default: 0
              params["nac_method"] - (0, 1, 2) - the expectations about what format to read
                     0 - non-SOC, non-polarized
                     1 - non-SOC, spin-polarized
                     2 - SOC, non-collinear
              params["minband"] - (int, starting from 1) - index of the lowest energy orbital to include 
                                  in the active space
              params["maxband"] - (int, starting from 1) - index of the highest energy orbital to include 
                                  in the active space
  
    Returns: The function returs lists containing: 
             e - energies
             coeff - MOs the plane wave coefficients,
             grid - the grid point vectors

             The number of elements in each list is determined by the number of k points
             Note that, for spin-polarized calculations, the number of k-points is always twice
             that of the non-spin-polarized or non-collinear k-points
    """  

    tim = Timer()
    tim.start()

    wd = get_value(params,"wd","wd","s")
    rd = get_value(params,"rd",os.getcwd()+"../../res","s") # of where the files will be printed out
    prefix = get_value(params,"prefix","x0.export","s")
    is_wfc = get_value(params,"read_wfc",1,"i")
    is_grd = get_value(params,"read_grid",1,"i")
    verb0 = get_value(params,"verb0",0,"i")
    verb1 = get_value(params,"verb1",0,"i")
    verb2 = get_value(params,"verb2",0,"i")
    nac_method = get_value(params,"nac_method",0,"i") 
    minband = get_value(params,"minband",1,"i")
    maxband = get_value(params,"maxband",2,"i")

    print "printing prefix:  ", prefix

    act_space = []    
    if(nac_method==0 or nac_method==1):
        file0 = "%s/%s/index.xml" % (wd, prefix)
        act_space = range(minband, maxband+1)              # min = 1, max = 2 => range(1,3) = [1,2]

        print "Reading index from file ", file0
        info, e = QE_methods.read_qe_index(file0, act_space, 1)

    elif nac_method==2:
        file1 = "%s/%s/index.xml" % (wd, prefix)
        act_space = range(2*minband-1, 2*(maxband+1)-1 )   # min =1, max = 2 => range(1,5) = [1,2,3,4]

        print "Reading index from file ", file1
        info, e = QE_methods.read_qe_index(file1, act_space, verb0)

    coeff = []
    grid = []

    for ik in xrange(info["nk"]):
        print "Handling the k-point %i with coordinates: %8.5f %8.5f %8.5f " \
         % (ik, info["k"][ik].x, info["k"][ik].y, info["k"][ik].z)

        if is_wfc==1:

            file2 = "%s/%s/wfc.%i" % (wd, prefix, ik+1)
            print "Reading the wfc from file ",file2
            coeff.append( QE_methods.read_qe_wfc(file2, act_space, verb1))   # CMATRIX(npw x len(act_space))

            ## Extract binary wfc file(s) ##
            #rd_wfc = rd + "/wf_storage_"+str(params["curr_index"])
            #pref1  = prefix[0:5]
            #pref2  = prefix[6:8]
            #os.system( "mkdir %s" % (rd_wfc) )
            #os.system( "cp %s/%s/%s.wfc* %s/" % (wd, pref1, pref2, rd_wfc) )

        if is_grd==1:
            file3 = "%s/%s/grid.%i" % (wd, prefix, ik+1)
            print "Reading the grid from file ", file3
            grid.append( QE_methods.read_qe_wfc_grid(file3 , verb2) )

    print "The time to read index, wavefunctions, and grid about your QE calculations = ", tim.stop();  

    return info, e, coeff, grid




def read_wfc_grid(params, info0, info1):
    """
    Read the coefficients and energies for the multi k-points cases, 
    even if some cases require gamma only

    \param[in] params A dictionary containing important simulation parameters
    \param[in] info0 Name of the temporary directory where data will be stored 
              for the case without SOC 
    \param[in] info1 Name of the temporary directory where data will be stored 
              for the case with SOC 
  
    Returns: Dictionaries that contain the data for the eigenvalues, g-vectors, and pw coefficients 
             for the current and next timesteps.    
    """

    tim = Timer()
    tim.start()


    nac_method = get_value(params,"nac_method",0,"i")  # choose what method for NAC calculations to use: 0 -standard, 1-corrected
    orthogonalize = get_value(params,"orthogonalize","0","i")

    """
    Here, adi - refers to spin-adiabatic (2-component spinor functions)
    Here, dia - refers to spin-diabatic (regular functions)
     
    """
    res_curr = {"Coeff_dia": None, "E_dia": None, "Coeff_adi": None, "E_adi":None, "grid": None}
    res_next = {"Coeff_dia": None, "E_dia": None, "Coeff_adi": None, "E_adi":None, "grid": None}

                
    if nac_method == 0 or nac_method == 1 or nac_method == 3:

        #====== Current electronic structure ========
        params["prefix"] = "curr0/x0.export"
        info_curr, e_curr, coeff_curr, grid_curr = read_all(params)

        if orthogonalize==1:
            print "Do internal orbital orthogonalization"
            coeff_curr[0] = orthogonalize_orbitals(coeff_curr[0])

            id1 = CMATRIX(coeff_curr[0].num_of_cols, coeff_curr[0].num_of_cols)
            id1.identity()
            if abs( (coeff_curr[0].H() * coeff_curr[0] - id1).max_elt() ) > 1e-5:
                print "Error\n"
                sys.exit(0)


        C_dia_curr, E_dia_curr = post_process(coeff_curr, e_curr, 0)
        res_curr["Coeff_dia"] = C_dia_curr
        res_curr["E_dia"] = E_dia_curr
        res_curr["grid"] = grid_curr

        #====== Next electronic structure ===========
        params["prefix"] = "next0/x0.export"
        info_next, e_next, coeff_next, grid_next = read_all(params)

        if orthogonalize==1:
            print "Do internal orbital orthogonalization"
            coeff_next[0] = orthogonalize_orbitals(coeff_next[0])

        C_dia_next, E_dia_next = post_process(coeff_next, e_next, 0)
        res_next["Coeff_dia"] = C_dia_next
        res_next["E_dia"] = E_dia_next
        res_next["grid"] = grid_next


    if nac_method == 2 or nac_method == 3:

        #====== Current electron electructure =======
        params["prefix"] = "curr1/x1.export"
        info_curr, e_curr, coeff_curr, grid_curr = read_all(params)

        if orthogonalize==1:
            print "Do internal orbital orthogonalization"
            coeff_curr[0] = orthogonalize_orbitals(coeff_curr[0])

            id1 = CMATRIX(coeff_curr[0].num_of_cols, coeff_curr[0].num_of_cols)
            id1.identity()
            if abs( (coeff_curr[0].H() * coeff_curr[0] - id1).max_elt() ) > 1e-5:
                print "Error\n"
                sys.exit(0)

        C_adi_curr, E_adi_curr = post_process(coeff_curr, e_curr, 1)
        res_curr["Coeff_adi"] = C_adi_curr
        res_curr["E_adi"] = E_adi_curr
        res_curr["grid"] = grid_curr

       
        #====== Next electronic structure ===========
        params["prefix"] = "next1/x1.export"
        info_next, e_next, coeff_next, grid_next = read_all(params)

        if orthogonalize==1:
            print "Do internal orbital orthogonalization"
            coeff_next[0] = orthogonalize_orbitals(coeff_next[0])

        C_adi_next, E_adi_next = post_process(coeff_next, e_next, 1)
        res_next["Coeff_adi"] = C_adi_next
        res_next["E_adi"] = E_adi_next
        res_next["grid"] = grid_next


    print "Time to read index, wfc, and wfc grids = ", tim.stop();

    return res_curr, res_next




def run(params):
    """
    This function is the main driver for workflows/pyxaid2

    \param[in] params A dictionary containing important simulation parameters
    """

    print "Starting trajectory.run"

    tim = Timer()
    # Parameters meaning
    # pp_type - pseudopotential type: US - ultra-soft, NC - norm-conserving, PAW - projector-augmented waves
    # wd - working directory, where all output (working) files will be written
    # rd - results directory, where all final results (energy, NAC, H', etc.) will be written by default it will be set to wd    

    # Now try to get parameters from the input
    critical_params = [ "wd"  ]
    default_params = { "BATCH_SYSTEM":"srun", "NP":1, "EXE":"", "EXE_EXPORT":"", "EXE_CONVERT":"",
                       "start_indx":0, "stop_indx":1, "dt":41.34145, "pp_type":"NC", 
                       "rd":params["wd"], "minband":1, "maxband":2, "minband_soc":1, "maxband_soc":2, 
                       "nac_method":0, "prefix0":"x0.scf", "prefix1":"x1.scf", "compute_Hprime":0 }
    comn.check_input(params, default_params, critical_params)


    BATCH_SYSTEM = params["BATCH_SYSTEM"]
    NP = params["NP"]
    EXE = params["EXE"]
    EXE_EXPORT = params["EXE_EXPORT"]
    EXE_CONVERT = params["EXE_CONVERT"]
    start_indx = int( params["start_indx"] )
    stop_indx = int( params["stop_indx"] )
    dt = params["dt"]
    pp_type = params["pp_type"]  
    wd = params["wd"]
    rd = params["rd"]
    minband = params["minband"]
    maxband = params["maxband"]
    minband_soc = params["minband_soc"]
    maxband_soc = params["maxband_soc"]
    nac_method = params["nac_method"]
    prefix0 = params["prefix0"]
    prefix1 = params["prefix1"]
    compute_Hprime = params["compute_Hprime"]


    # Sanity/Convention check
    if(minband<=0): 
        print "Error: minband should be >0, current value of minband = ",minband
        sys.exit(0)
    if(minband>maxband):
        print "Error: minband must be smaller or equal to maxband. Current values: minband = ",minband," maxband = ",maxband
        sys.exit(0)
    
    if nac_method == 0:
        print "non-relativistic, non-spin-polarized calculations \n"
    elif nac_method == 1:
        print "non-relativistic, spin-polarized calculations \n"
    elif nac_method == 2:
        print "fully relativistic, non-collinear calculations \n"
    elif nac_method == 3:
        print "same as 2 + 1, followed by the projections  \n"
    else:
        print "Error: nac_method must be one of the values in [0,1,2,3]  \n"
        sys.exit(0)


    # Use this for nspin = 1 or 2
    act_sp1 = range(minband, maxband+1)     # min = 1, max = 2 => range(1,3) = [1,2]

    # Use this for nspin = 4
    act_sp2 = range(2*minband_soc-1, 2*(maxband_soc+1)-1 ) # min =1, max = 2 => range(1,5) = [1,2,3,4]



    # Initialize variables
    curr_index = start_indx - 1
    t = start_indx

    print "In trajectory.run: current working directory for python: ",os.getcwd()
    print "In trajectory.run: current working directory for sh:",os.system("echo $(pwd)")

    # Create the working directory where all output files will be written
    # The results directory should already exist
    os.system("mkdir %s" % wd)  
    while t<=stop_indx:
        print ">>>>>>>>>>>>>>>>>>>>  t= ", t, " <<<<<<<<<<<<<<<<<<<<<"

        print "stop_indx", stop_indx
        print "start_indx", start_indx

        dirname = ""
        if t==start_indx:
           print "Starting first point in this batch"
           dirname0, dirname1 = "curr0", "curr1"

        if t>start_indx:
           print "Continuing with other points in this batch"
           dirname0, dirname1 = "next0", "next1"


        # Run the QE calculations       
        run_qe(params, t, dirname0, dirname1)

        
        if curr_index>=start_indx:

            # Update curr_index in params - this may be useful
            params["curr_index"] = curr_index

            # First see wether the calculation is what we wanted
            info0, all_e_dum0, info1, all_e_dum1 = read_info(params)

            # Read present and next electronic structure (es) in a formatted way
            es_curr, es_next = read_wfc_grid(params, info0, info1)

            # Finally,  using the current and the next wavefunctions to compute the properties of interest
            # non-relativistic, non-spin-polarized case
            if nac_method == 0 or nac_method == 1 or nac_method == 3:

                if info0["nspin"]==1 or info0["nspin"]==2:  # non SOC case

                    if info0["nk"]==1 or info0["nk"] == 2: # Only one k-point

                        H_nosoc = compute_properties.compute_properties_dia_gamma(params, es_curr, es_next, curr_index)

                    else: 
                        compute_properties.compute_properties_general(params, es_curr, es_next, curr_index)
           
                    if compute_Hprime == 1:
                        compute_hprime.compute_hprime_dia(es_curr, info0, "%s/0_Hprime_%d" % (rd, curr_index) )
               

            if nac_method == 2 or nac_method == 3:

                if info1["nk"]==1: # Only one k-point
               
                    H_soc = compute_properties.compute_properties_adi_gamma(params, es_curr, es_next, curr_index)

                else:
                    print "Multiple k-points scheme with SOC is not yet implemented"
                    sys.exit(0)

            
            # spin-polarized case
            if  nac_method == 3:
                # check whether the adiabatic and diabatic basis have the same number of plane waves
                # the reason why I used the read_qe_wfc_info is because I will need the ngw 
                # to check the consistency 
                # But the read_qe_index does not read it, so in order to avoid the changes in the Libra code, 
                # I use the read_qe_wfc_info.
                info_wfc0 = QE_methods.read_qe_wfc_info("%s/curr0/x0.export/wfc.1" % wd,0)
                info_wfc1 = QE_methods.read_qe_wfc_info("%s/curr1/x1.export/wfc.1" % wd,0)

                if info_wfc0["ngw"] != info_wfc1["ngw"]:
                    print "Error: the number of plane waves if diabatic and adiabatic functions should be equal"
                    sys.exit(0)


                params = {"do_orth": 0, "root_directory": rd, "curr_index": curr_index, "print_overlaps": 1, "dt": dt}
                compute_ovlps(coeff_curr0, coeff_next0, coeff_curr1, coeff_next1, e_curr0, e_next0, e_curr1, e_next1, params)

            """
            if nac_method == 0 or nac_method == 1 or nac_method == 3:
            
                H_nosoc.real().show_matrix("%s/0_Ham_%d_re" % (rd, curr_index) )
                H_nosoc.imag().show_matrix("%s/0_Ham_%d_im" % (rd, curr_index) )
            
            if nac_method == 2 or nac_method == 3:
                H_soc.real().show_matrix("%s/0_Ham_soc_%d_re" % (rd, curr_index) )
                H_soc.imag().show_matrix("%s/0_Ham_soc_%d_im" % (rd, curr_index) )
            """

            #-----------------------------------------------------------------

            # Remove current run, make next run to be the current one
            if nac_method == 0 or nac_method == 1 or nac_method==3:
                os.system("rm -rf %s/curr0" % wd )
                os.system("mv %s/next0 %s/curr0" % (wd, wd) )
            
            if nac_method==2 or nac_method==3:
                os.system("rm -rf %s/curr1" % wd )
                os.system("mv %s/next1 %s/curr1" % (wd, wd) )
            
            print "old files deleted, new have become old"
             
            
# ACHTUNG!!! Restoring wfc makes some complications, so we might need to destroy wfc objects
# after each round of operations and create new objects from the beginning - thia may be safer!

        curr_index = curr_index + 1
               
        print "End of step t=", t
        t = t + 1






