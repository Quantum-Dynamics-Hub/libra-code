#***********************************************************
# * Copyright (C) 2017-2019 Brendan A. Smith, Wei Li, and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

"""
.. module:: step2
   :platform: Unix, Windows
   :synopsis: This module implements functions for building molecular structures

.. moduleauthor:: Brendan A. Smith, Wei Li, and Alexey V. Akimov

"""


import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

from utils import *
import libra_py.common_utils as comn
from libra_py import units
import compute_properties
import compute_hprime


def run_qe(params, t, dirname0, dirname1):
    """

    This function runs necessary QE calculations as defined by the "params" dictionary

    Args:
        params ( dictionary ): A dictionary containing important simulation parameters

            * **params["BATCH_SYSTEM"]** ( string ): the name of the job submission command
                use "srun" if run calculations on SLURM system or "mpirun" if run on PBS system
                [default: "srun"]
            * **params["NP"]** ( int ): the number of nodes on which execute calculations
                [default: 1]
            * **params["EXE"]** ( string ): the name of the program to be executed. This may be 
                the absolute path to the QE (pw.x) binary
            * **params["EXE_EXPORT"]** ( string ): the name of the program that converts the binary files
                with the QE wavefunctions to the text format (pw_export.x). The name includes the 
                absolute path to the binary
            * **params["prefix0"]** ( string ): the name of scf template input file - it should 
                contain all the parameters controlling the computational methodology for QE.
                If the file is called "x0.scf.in", use "x0.scf" as the value of the "prefix0"
                [default: "x0.scf"]
            * **params["prefix1"]** ( string ): the name of scf template input file - it should 
                contain all the parameters controlling the computational methodology for QE.
                Presently is used for SOC-enabled calculations, whereas the "prefix0" defines the
                no-SOC calculations. If the file is called "x1.scf.in", use "x1.scf" as the value
                of the "prefix1" [default: "x1.scf"]
            * **params["nac_method"]** ( int ): selects the type of calculations to perform:
 
                - 0 : non-spin-polarized calculations (needs only "prefix0")
                - 1 : spin-polarized calculations (needs only "prefix0")
                - 2 : non-collinear calculation (SOC) only (needs only "prefix1")
                - 3 : spin-polarized and non-collinear calculation (SOC) (needs both "prefix0" and "prefix1")

                [default: 0]

            * **params["wd"]** ( string ): the name of a "working directory (can be removed once the calculatons
                are done)" that will be created during this function execution.

        t ( int ): the current time step
        dirname0 ( string ): Name of the temporary directory where data will be stored 
            for the case without SOC 
        dirname1 ( string ): Name of the temporary directory where data will be stored 
            for the case with SOC 

    """

    tim = Timer()
    tim.start()


    # Now try to get parameters from the input
    critical_params = [ "EXE", "EXE_EXPORT" ] 
    default_params = { "BATCH_SYSTEM":"srun", "NP":1, "prefix0":"x0.scf", "prefix1":"x1.scf", "nac_method":0, "wd":"wd"  }
    comn.check_input(params, default_params, critical_params)


    BATCH_SYSTEM = params["BATCH_SYSTEM"]
    NP = params["NP"]
    EXE = params["EXE"]
    EXE_EXPORT = params["EXE_EXPORT"]
    prefix0 = params["prefix0"]
    prefix1 = params["prefix1"]
    nac_method = params["nac_method"]
    wd = params["wd"]

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

    This fucntions reads the output from QE calculations, and stores the output
    information in dictionaries

    Args:
        params ( dictionary ): Calculation control parameters
        
            * **params["nac_method"]** ( int ): selects the type of output to analyze:

                - 0 : non-spin-polarized calculations 
                - 1 : spin-polarized calculations
                - 2 : non-collinear calculation (SOC) only 
                - 3 : spin-polarized and non-collinear calculation (SOC)
  
    Returns: 
    
        tuple: ( info0, all_e_dum0, info1, all_e_dum1 ): 

            info0 ( dictionary ): QE calculations info for the spin-diabatic calculations

            all_e_dum0 ( list of CMATRIX(norb, norb) objects ): (eigen)energies for all the k-points for 
                the spin-diabatic calculations

            info1 ( dictionary ): QE calculations info for the non-collinear (spin-adiabatic) calculations

            all_e_dum1 ( list of CMATRIX(norb, norb) objects ): (eigen)energies for all the k-points for 
                the non-collinear (spin-adiabatic) calculations

            ..seealso:: ```QE_methods.read_qe_index```

    """

    tim = Timer()
    tim.start()

    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "nac_method":0, "wd":"wd" }
    comn.check_input(params, default_params, critical_params)

    nac_method = params["nac_method"] 
    wd0 = params["wd"] 



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

    Args:

    params ( dictionary ): Parameters controlling the simulation parameters

        * **params["prefix"]** ( string ): the location of the folder containing index.xml, wfc.*, and grid.* files [default: x0.export ]
        * **params["read_wfc"]** ( 0 or 1 ): whether or not to read the wfc coefficients. [ default: 1 ]
        * **params["read_grid"]** ( 0 or 1 ): whether or not to read the grid informations. [ default: 1 ]
        * **params["verb0"]** ( 0 or 1 ): turn off/on the extra printout while reading index.xml. [ default: 0 ]
        * **params["verb1"]** ( 0 or 1 ): turn off/on the extra printout while reading wfc.*. [ default: 0 ]
        * **params["verb2"]** ( 0 or 1 ): turn off/on the extra printout while reading grid.*. [ default: 0 ]
        * **params["nac_method"]** ( 0, 1, 2 ): the expectations about what format to read:

            - 0 - non-SOC, non-polarized
            - 1 - non-SOC, spin-polarized
            - 2 - SOC, non-collinear

        * **params["minband"]** ( int ): index of the lowest energy orbital to include
            in the active space, counting starts from 1 [ default: 1]

        * **params["maxband"]** ( int ): index of the highest energy orbital to include 
            in the active space, counting starts from 1 [ defaults: 2]
  
    Returns: 
        tuple: ( info, e, coeff, grid ), where 

            * info ( dictionary ): general descritor info ..seealso::```QE_methods.read_qe_index```
            * e ( list of CMATRIX(norbs, norbs) ): band energies for each k-pints  ..seealso::```QE_methods.read_qe_index```
            * coeff ( list of CMATRIX(npw, len(act_space)) objects ): such the 
                coeff[k] are the MOs in the plane wave basis for the k-point k
            * grid ( list of VECTOR objects ): the grid point vectors [ units: tpiba ]

            The number of elements in each list is determined by the number of k points
            Note that, for spin-polarized calculations, the number of k-points is always twice
            that of the non-spin-polarized or non-collinear k-points
    """  

    tim = Timer()
    tim.start()

    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "nac_method":0, "wd":"wd" , "rd":os.getcwd()+"../../res",  "prefix":"x0.export" 
    "read_wfc":1, "read_grid":1, "verb0":0, "verb1":0, "verb2":0, "minband":1, "maxband":2  }
    comn.check_input(params, default_params, critical_params)

    wd = params["wd"]
    rd = params["rd"]
    prefix = params["prefix"]
    is_wfc = params["read_wfc"]
    is_grd = params["read_grid"]
    verb0 = params["verb0"]
    verb1 = params["verb1"]
    verb2 = params["verb2"]
    nac_method = params["nac_method"]
    minband = params["minband"]
    maxband = params["maxband"]


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




def read_wfc_grid(params):
    """

    Read the coefficients and energies for the multi k-points cases, 
    even if some cases require gamma only

    Args:
        params ( dictionary ): The control parameters, which may contain:

            * **param["nac_method"]** ( 0, 1, 2, 3 ): the expectations about what format to read:

                - 0 : non-spin-polarized calculations [ default ]
                - 1 : spin-polarized calculations
                - 2 : non-collinear calculation (SOC) only 
                - 3 : spin-polarized and non-collinear calculation (SOC)
         
  
    Returns: 
        tuple: ( res_curr, res_next ), Here _curr, refers to the current timestep properties,
            and the _next refers to the consecutive timestep properties. Each element of the 
            output is a dictionary with the following elements:

            * **res_curr["Coeff_dia"]** ( list of CMATRIX(npw, len(act_space)) objects ) : the 
                wavefunction coefficients in the planewave basis for the spin-diabatic wavefunctions, such
                that res_curr["Coeff_dia"][k] is a matrix for the k-point with index k.
                Only for nac_method == 0, 1, and 3. 

            * **res_curr["E_dia"]** ( list of MATRIX(len(act_space), len(act_space)) objects ) : the MO
                energies for the spin-diabatic wavefunctions, such
                that res_curr["E_dia"][k] is a matrix for the k-point with index k.
                Only for nac_method == 0, 1, and 3. 

            * **res_curr["Coeff_adi"]** ( list of CMATRIX(npw, len(act_space)) objects ) : the 
                wavefunction coefficients in the planewave basis for the spin-adiabatic wavefunctions, such
                that res_curr["Coeff_adi"][k] is a matrix for the k-point with index k.
                Only for nac_method == 2 and 3. 

            * **res_curr["E_adi"]** ( list of MATRIX(len(act_space), len(act_space)) objects ) : the MO
                energies for the spin-adiabatic wavefunctions, such
                that res_curr["E_adi"][k] is a matrix for the k-point with index k.
                Only for nac_method == 2 and 3. 

            * **res_curr["grid"]** ( list of VECTOR objects ): the grid point vectors [ units: tpiba ]

    """

    tim = Timer()
    tim.start()


    # Now try to get parameters from the input
    critical_params = [ ] 
    default_params = { "nac_method":0, "orthogonalize":0 }
    comn.check_input(params, default_params, critical_params)

    nac_method = params["nac_method"]
    orthogonalize = params["orthogonalize"]



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
    """This function is the main driver for the NAC calculations withing the workflows/nbra

    Args:
        params ( dictionary ): Simulation control parameters

        * **params["start_indx"]** ( int ): index of the starting datapoint in the trajectory to compute
            the NACs [default: 0]
        * **params["stop_indx"]** ( int ): index of the final datapoint in the trajectory to compute
            the NACs [default: 1]
        * **params["dt"]** ( double ): time step between two adjacent datapoint a the trajectory [units: a.u.; default: 41.0]
        * **params["wd"]** ( string ): the name of a "working directory (can be removed once the calculatons
            are done)" that will be created during this function execution - this is where all output 
            (working) files will be written [ default: "wd" ]
        * **params["nac_method"]** ( int ): selects the type of output to analyze:

            - 0 : non-spin-polarized calculations [default]
            - 1 : spin-polarized calculations
            - 2 : non-collinear calculation (SOC) only 
            - 3 : spin-polarized and non-collinear calculation (SOC)

        * **params["minband"]** ( int ): index of the lowest energy orbital to include
            in the active space in the non-SOC (spin-diabatic) calculations, 
            counting starts from 1. Used for nac_method == 0, 1, 3. [ default: 1]
        * **params["maxband"]** ( int ): index of the highest energy orbital to include 
            in the active space  in the non-SOC (spin-diabatic) calculations, 
            counting starts from 1. Used for nac_method == 0, 1, 3. [ defaults: 2]
        * **params["minband_soc"]** ( int ): index of the lowest energy orbital pair to include
            in the active space in the SOC (spin-adiabatic) calculations, 
            counting starts from 1. Used for nac_method == 2 and 3. [ default: 1]
        * **params["maxband_soc"]** ( int ): index of the highest energy orbital pair to include 
            in the active space  in the SOC (spin-adiabatic) calculations, 
            counting starts from 1. Used for nac_method == 2 and 3. [ defaults: 2]
        * **params["compute_Hprime"]** ( Boolean ): the flag to compute the <i|H'|j> matrices [ default: False]
        * **params["verbosity"]** ( int ): the verbosity level regarding the execution of the current function [default: 0]

    Returns:
        None: but generates the files with couplings, energies and transition density matrices
 
    """

    # Now try to get parameters from the input
    critical_params = [ ]
    default_params = { "start_indx":0, "stop_indx":1, "dt":1.0*units.fs2au, 
                       "wd":"wd",
                       "nac_method":0,
                       "minband":1, "maxband":2, 
                       "minband_soc":1, "maxband_soc":2, 
                       "compute_Hprime":False, 
                       "verbosity":0
                     }
    comn.check_input(params, default_params, critical_params)


    start_indx =  params["start_indx"] 
    stop_indx = params["stop_indx"] 
    dt = params["dt"]
    wd = params["wd"]
    nac_method = params["nac_method"]
    minband = params["minband"]
    maxband = params["maxband"]
    minband_soc = params["minband_soc"]
    maxband_soc = params["maxband_soc"]
    compute_Hprime = params["compute_Hprime"]
    verbosity = params["verbosity"]

    if verbosity>0:
        print "Starting trajectory.run"
    tim = Timer()

    # Sanity/Convention check
    if(minband<=0): 
        print "Error: minband should be >0, current value of minband = ",minband
        sys.exit(0)
    if(minband>maxband):
        print "Error: minband must be smaller or equal to maxband. Current values: minband = ",minband," maxband = ",maxband
        sys.exit(0)
        
    if nac_method == 0:        
        if verbosity>0:
            print "non-relativistic, non-spin-polarized calculations \n"
    elif nac_method == 1:
        if verbosity>0:
            print "non-relativistic, spin-polarized calculations \n"
    elif nac_method == 2:
        if verbosity>0:
            print "fully relativistic, non-collinear calculations \n"
    elif nac_method == 3:
        if verbosity>0:
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

    if verbosity>0:
        print "In trajectory.run: current working directory for python: ",os.getcwd()
        print "In trajectory.run: current working directory for sh:",os.system("echo $(pwd)")

    # Create the working directory where all output files will be written
    # The results directory should already exist
    os.system("mkdir %s" % wd)  
    while t<=stop_indx:

        if verbosity>1:
            print ">>>>>>>>>>>>>>>>>>>>  t= ", t, " <<<<<<<<<<<<<<<<<<<<<"

            print "stop_indx", stop_indx
            print "start_indx", start_indx

        dirname = ""
        if t==start_indx:
           if verbosity>1:
               print "Starting first point in this batch"
           dirname0, dirname1 = "curr0", "curr1"

        if t>start_indx:
           if verbosity>1:
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
           
                    if compute_Hprime == True:
                        compute_hprime.compute_hprime_dia(es_curr, info0, "%s/0_Hprime_%d" % (wd, curr_index) )
               

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


                params1 = {"do_orth": 0, "root_directory": rd, "curr_index": curr_index, "print_overlaps": 1, "dt": dt}
                compute_ovlps(coeff_curr0, coeff_next0, coeff_curr1, coeff_next1, e_curr0, e_next0, e_curr1, e_next1, params1)

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
            
            if verbosity>1:
                print "old files deleted, new have become old"
             
            
        # ACHTUNG!!! Restoring wfc makes some complications, so we might need to destroy wfc objects
        # after each round of operations and create new objects from the beginning - thia may be safer!

        curr_index = curr_index + 1
               
        if verbosity>1:
            print "End of step t=", t

        t = t + 1






