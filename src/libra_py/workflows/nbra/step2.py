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
   :synopsis: This module implements functions for doing NAC calculaitons in the 1-electron
       basis of KS orbitals, withe the QE package

.. moduleauthor:: Brendan A. Smith, Wei Li, and Alexey V. Akimov

"""

import os
import sys
# The sparse matrix library for storing the MO overlaps
# with a high number of states
import numpy as np
import scipy.sparse
import time
# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import libra_py.packages.qe.methods as QE_methods
import libra_py.packages.qe.utils as QE_utils
import libra_py.packages.cp2k.methods as CP2K_methods
from libra_py import cube_file_methods
from libra_py import molden_methods
from libra_py import data_conv
from libra_py import units

#import libra_py.common_utils as comn
import util.libutil as comn

from . import compute_properties
from . import compute_hprime



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
            are done)" that will be created during this function execution - this is where the temporary files
            are written to 
        * **params["rd"]** ( string ): the name of a "results directory" - this is where all the output files 
            files will be written to [ default: "res" ]
        * **params["nac_method"]** ( int ): selects the type of output to analyze:

            - 0 : non-spin-polarized calculations [default]
            - 1 : spin-polarized calculations
            - 2 : non-collinear calculation (SOC) only 
            - 3 : spin-polarized and non-collinear calculation (SOC)

        * **params["compute_Hprime"]** ( Boolean ): the flag to compute the <i|H'|j> matrices [ default: False]
        * **params["verbosity"]** ( int ): the verbosity level regarding the execution of the current function [default: 0]


        ..seealso:: parameters of ::funct::```run_qe``` function
            * **params["BATCH_SYSTEM"]**
            * **params["NP"]**
            * **params["EXE"]**
            * **params["EXE_EXPORT"]**
            * **params["prefix0"]**
            * **params["prefix1"]**

        ..seealso:: parameters of ::funct::```read_all``` function

            * **params["read_wfc"]** 
            * **params["read_grid"]**
            * **params["verb0"]** 
            * **params["verb1"]**
            * **params["verb2"]**

        ..seealso:: parameters of ::funct::```read_wfc_grid``` function
            * **params["maxband"]**
            * **params["minband"]**
            * **params["maxband_soc"]**
            * **params["minband_soc"]**



    Returns:
        None: but generates the files with couplings, energies and transition density matrices
            The matrices are of  2N x 2N dimensions, where N = maxband - minband + 1
 
    """

    # Now try to get parameters from the input
    critical_params = [ "rd" ]
    default_params = { "start_indx":0, "stop_indx":1, "dt":1.0*units.fs2au, 
                       "wd":"wd",
                       "nac_method":0,
                       "pdos_flg":0,
                       "compute_Hprime":False, 
                       "verbosity":0
                     }
    comn.check_input(params, default_params, critical_params)


    rd = params["rd"]

    start_indx = int(params["start_indx"]) 
    stop_indx  = int(params["stop_indx"]) 
    dt = params["dt"]
    wd = params["wd"]
    nac_method = params["nac_method"]
    compute_Hprime = params["compute_Hprime"]
    verbosity = params["verbosity"]    
    pdos_flg = params["pdos_flg"]

    if verbosity>0:
        print( "Starting trajectory.run")
    tim = Timer()

    # Sanity/Convention check
        
    if nac_method == 0:        
        if verbosity>0:
            print( "non-relativistic, non-spin-polarized calculations \n")
    elif nac_method == 1:
        if verbosity>0:
            print( "non-relativistic, spin-polarized calculations \n")
    elif nac_method == 2:
        if verbosity>0:
            print( "fully relativistic, non-collinear calculations \n")
    elif nac_method == 3:
        if verbosity>0:
            print( "same as 2 + 1, followed by the projections  \n")
    else:
        print( "Error: nac_method must be one of the values in [0,1,2,3]  \n")
        sys.exit(0)



    # Initialize variables
    curr_index = start_indx - 1
    t = start_indx

    if verbosity>0:
        print( "In trajectory.run: current working directory for python: ",os.getcwd())
        print( "In trajectory.run: current working directory for sh:",os.system("echo $(pwd)"))

    # Create the working directory where all output files will be written
    # The results directory should already exist
    os.system("mkdir %s" % wd)  
    while t<=stop_indx:

        if verbosity>0:
            print( ">>>>>>>>>>>>>>>>>>>>  t= ", t, " <<<<<<<<<<<<<<<<<<<<<")

            print( "stop_indx", stop_indx)
            print( "start_indx", start_indx)

        dirname = ""
        if t==start_indx:
           if verbosity>0:
               print( "Starting first point in this batch")
           dirname0, dirname1 = "curr0", "curr1"

        if t>start_indx:
           if verbosity>0:
               print( "Continuing with other points in this batch")
           dirname0, dirname1 = "next0", "next1"


        # Run the QE calculations       
        QE_methods.run_qe(params, t, dirname0, dirname1)

        
        if curr_index>=start_indx:

            # Update curr_index in params - this may be useful
            params["curr_index"] = curr_index

            # First see wether the calculation is what we wanted
            info0, all_e_dum0, info1, all_e_dum1 = QE_methods.read_info(params)

            # Read present and next electronic structure (es) in a formatted way
            es_curr, es_next = QE_methods.read_wfc_grid(params) 

            # Finally,  using the current and the next wavefunctions to compute the properties of interest

            # For the non-relativistic and either spin or non-spin-polarized case

            if nac_method == 0 or nac_method == 1 or nac_method == 3:
                if ( info0["nspin"]==1 and info0["nk"]==1 ) or ( info0["nspin"]==2 and info0["nk"] == 2) : # Non-SOC case

                    # Only one k-point. For the spin-polarized case, the beta orbtials are as if they were computed using 2 K-points
                    # (nk = 2). However, they are not. This is just the QE format for expressing alpha and beta orbitals at a single K-point 
                    #if info0["nk"]==1 or info0["nk"] == 2:  

                    if verbosity>0:
                        print( "Computing various properies of the spin-diabatic (non-relativistic) KS orbitals using a single K-point")

                    # Compute various properties of the spin-diabatic (non-relativistic) KS orbitals at a single K-point 
                    compute_properties.compute_properties_onekpt(params, es_curr, es_next, curr_index)

                else: 

                    # Compute various properies of the spin-diabatic (non-relativistic) KS orbitals using multiple K-point
                    if verbosity>0:
                        print( "Computing various properies of the spin-diabatic (non-relativistic) KS orbitals using multiple K-points")
                        print( "Warning: This capabilitiy is under development and not fully tested")

                    print ("This capability is temporarily disabled. Please use only a single K-point for now")
                    print ("Exiting now")
                    sys.exit(0)
                    #compute_properties.compute_properties_general(params, es_curr, es_next, curr_index)
           
                if compute_Hprime == True:

                    # Compute the transition dipole moment along the nuclear trajectory
                    if verbosity>0:
                        print( "Computing the transition dipole moment along the nuclear trajectory")
                        print( "Warning: This capabilitiy is under development and not fully tested")

                    # C++ implementation
                    compute_hprime.compute_hprime_dia(es_curr, info0, "%s/0_Hprime_%d" % (rd, curr_index) )

                    # Python implementation
                    #compute_hprime.hprime_py(es_curr, info0, "%s/0_Hprime_%d" % (rd, curr_index) )
               

            # For the relativistic case
            if nac_method == 2 or nac_method == 3:

                # Only one k-point
                if info1["nk"]==1:

                    if verbosity>0:
                        print( "Computing various properies of the spin-adiabatic (relativistic) KS orbitals using a single K-point")

                    # Compute various properties of the spin-adiabatic (relativistic) KS orbitals at a single K-point              
                    compute_properties.compute_properties_onekpt(params, es_curr, es_next, curr_index)
                    #compute_properties.compute_properties_adi_gamma(params, es_curr, es_next, curr_index)

                else:
                    print( "Multiple k-points scheme with SOC is not yet implemented")
                    sys.exit(0)

            
            # Checking if spin-adiabtic and spin-diabatic cases were using the same sized P\plane wave basis
            if  nac_method == 3:

                print( "nac_method == 3: Entering check for whether the adiabatic and diabatic basis have the same number of plane waves")
             
                # Check whether the adiabatic and diabatic basis have the same number of plane waves
                # the reason why I used the read_qe_wfc_info is because I will need the ngw 
                # to check the consistency 
                # But the read_qe_index does not read it, so in order to avoid the changes in the Libra code, 
                # I use the read_qe_wfc_info.

                info_wfc0 = QE_methods.read_qe_wfc_info("%s/curr0/x0.export/wfc.1" % wd,0)
                info_wfc1 = QE_methods.read_qe_wfc_info("%s/curr1/x1.export/wfc.1" % wd,0)

                if info_wfc0["ngw"] != info_wfc1["ngw"]:
                    print( "Error: The number of plane waves of diabatic and adiabatic functions should be equal")
                    sys.exit(0)
                else:
                    print( "Pass: The number of plane waves of diabatic and adiabatic functions are equal")
 
                params1 = {"do_orth": 0, "root_directory": rd, "curr_index": curr_index, "print_overlaps": 1, "dt": dt}
                compute_ovlps(coeff_curr0, coeff_next0, coeff_curr1, coeff_next1, e_curr0, e_next0, e_curr1, e_next1, params1)

            #-----------------------------------------------------------------

            # Remove current run, make next run to be the current one
            if nac_method == 0 or nac_method == 1 or nac_method==3:
                os.system("rm -rf %s/curr0" % wd )
                os.system("mv %s/next0 %s/curr0" % (wd, wd) )
            
            if nac_method==2 or nac_method==3:
                os.system("rm -rf %s/curr1" % wd )
                os.system("mv %s/next1 %s/curr1" % (wd, wd) )
            
            if verbosity>0:
                print( "old files deleted, new have become old")
       
      
        # Move the index files to the results directory, just don't copy the very last one - it is a repetition
        #if t<stop_indx:
        if nac_method == 0 or nac_method == 1 or nac_method==3:
            # The input file template defined in params["prefix0"] should have   prefix = 'x0' !
            os.system("cp %s/curr0/x0.export/index.xml %s/x0_index_%i.xml" % (wd, rd, t))
            os.system("cp %s/curr0/x0.save/data-file-schema.xml %s/x0_data-file-schema_%i.xml" % (wd, rd, t))
            os.system("cp %s/curr0/%s.%d.out %s/%s.%d.out" % (wd, params["prefix0"], t, rd, params["prefix0"], t))

            if pdos_flg == 1:
                print ("Entering pdos_flag == 1: Outputting files necessary for computing pDOS calculations with projwfc.x")
                print ("Warning: Crashes at this point may result from QE versioning. This flag is currently compatable with QE v6.2.1")
                print ("Warning: This flag is tested only for non-relativisitic cases")

                os.system("cp %s/curr0/x0.save/charge-density.dat %s/x0_charge-density_%i.dat" % (wd, rd, t))
                os.system("cp %s/curr0/x0.save/wfcdw1.dat %s/x0_wfcdw1_%i.dat" % (wd, rd, t))
                os.system("cp %s/curr0/x0.save/wfcup1.dat %s/x0_wfcup1_%i.dat" % (wd, rd, t))

        if nac_method == 2 or nac_method == 3:
            # The input file template defined in params["prefix1"] should have   prefix = 'x1' !
            os.system("cp %s/curr1/x1.export/index.xml %s/x1_index_%i.xml" % (wd, rd, t))
            os.system("cp %s/curr1/x1.save/data-file-schema.xml %s/x1_data-file-schema_%i.xml" % (wd, rd, t))
            os.system("cp %s/curr1/%s.%d.out %s/%s.%d.out" % (wd, params["prefix1"], t, rd, params["prefix1"], t))
      
    
        # ACHTUNG!!! Restoring wfc makes some complications, so we might need to destroy wfc objects
        # after each round of operations and create new objects from the beginning - thia may be safer!

        curr_index = curr_index + 1
               
        if verbosity>1:
            print( "End of step t=", t)

        t = t + 1

def run_cp2k_libint_step2(params):
    """
    This function runs the step2 for computing the MO overlaps for both DFT and xTB calculations and saves them as sparse 
    matrices using scipy.sparse library. In order to load the saved sparse matrices you need to use
    scipy.sparse.load_npz command.

    Args:

        params (dictionary): A dictionary containing the following parameters:

                             res_dir (string): The directory for saving the MO overlaps.

                             all_logfiles (string): The directory to save all log files.

                             all_pdosfiles (string): The directory to save all pdos files.

                             istep (integer): The initial step.

                             fstep (integer): The final step.

                             lowest_orbital (integer): The lowest orbital considered in the calculations.

                             highest_orbital (integer): The highest orbital considered in the calcualtions.

                             is_spherical (bool): Flag for spherical coordinates.

                             isxTB (bool): Flag for xTB calculations.

                             isUKS (bool): Flag for unrestricted spin calculations.

                             remove_molden (bool): Flag for removing the molden files after computing the MO overlaps.

                             nprocs (integer): Number of processors.

                             cp2k_ot_input_template (string): The CP2K OT input template file name (just for the xTB).

                             cp2k_diag_input_template (string): The CP2K diagonalization input template file name. For
                                                                DFT calculations, this is the only parameter that we need.

                             trajectory_xyz_filename (string): The trajectory xyz file name.

                             cp2k_exe (string): The full path to CP2K executable.

                             mpi_executable (string): The mpi run executable for CP2K. It can be mpiexe, mpirun, or srun.

                             cube_visualization (bool): Flag for plotting the molecular orbitals using VMD.
 
                             vmd_input_template (string): The name of the VMD input template for visualizing the cube files.

                             states_to_plot (list): The list of states that need to be plot.

                             plot_phase_corrected (bool): Flag for plotting the molecular orbitals phase-corrected for the running job.

                             vmd_exe (string): The VMD executable.

                             tachyon_exe (string): The VMD Tachyon executable for rendering high-quality images.

                             x_pixels (integer): Number of pixels in the X direction of the image.

                             y_pixels (integer): Number of pixels in the Y direction of the image.

                             remove_cube (bool): Flag for removing the cube files after plotting the molecular orbitals.

    """

    # Now try to get parameters from the input
    critical_params = [ "cp2k_ot_input_template", "cp2k_diag_input_template", "trajectory_xyz_filename" ]
    default_params = {"res_dir": os.getcwd() + "/res", "all_logfiles": os.getcwd() + "/all_logfiles", "all_pdosfiles": os.getcwd() + "/all_pdosfiles", "all_images": os.getcwd() + "/all_images", "image_format": 'bmp', "istep": 0, "fstep": 2, "lowest_orbital": 1, "highest_orbital": 2, "is_spherical": True, "isXTB": False, "isUKS": False, "remove_molden": True, "nprocs": 2, "cp2k_exe": "cp2k.psmp", "mpi_executable": "mpirun", "cube_visualization": False, "vmd_input_template": "vmd.tcl", "states_to_plot": [1], "plot_phase_corrected": True, "vmd_exe": "vmd", "tachyon_exe": "tachyon_LINIXAMD64", "x_pixels": 1024, "y_pixels": 1024, "remove_cube": True, 'together_mode': False}
    comn.check_input(params, default_params, critical_params)


    # Making the required directories for gatherig all the information needed
    # including the overlap data and pdos files. The logfiles do not contain
    # specific data but we finally move them to all_logfiles
    if not os.path.exists(params['res_dir']):
        os.system(F"mkdir {params['res_dir']}")
    if not os.path.exists(params['all_logfiles']):
        os.system(F"mkdir {params['all_logfiles']}")
    if not os.path.exists(params['all_pdosfiles']):
        os.system(F"mkdir {params['all_pdosfiles']}")
    if not os.path.exists(params['all_images']):
        os.system(F"mkdir {params['all_images']}")

    # setting up the initial step and final step
    istep = params['istep']
    fstep = params['fstep']

    # spherical or cartesian GTOs
    is_spherical = params['is_spherical']

    # number of processors
    nprocs = params['nprocs']

    # CP2K executable
    cp2k_exe = params['cp2k_exe']

    # mpirun executable
    mpi_executable = params['mpi_executable']

    # Unrestricted spin calculations
    isUKS = params['isUKS']

    # Extended tight-binding calculations
    isxTB = params['isxTB']

    # Periodic calculation falg
    is_periodic = params['is_periodic']

    # Cube visualization
    cube_visualization = params['cube_visualization']

    if cube_visualization:
        # Read the tcl file data
        states_to_plot = params['states_to_plot']
        if isUKS:
            phase_factors_alpha = np.ones(( len(states_to_plot), 1 ))
            phase_factors_beta = np.ones(( len(states_to_plot), 1 ))
        else:
            phase_factors_alpha = np.ones(( len(states_to_plot), 1 ))

    # the counter for a job steps, this counter is needed for not reading
    # the data of a molden file twice
    counter = 0
    print('-----------------------Start-----------------------')
    for step in range(istep, fstep):
        # a timer for all the procedure
        t1_all = time.time()
        print('-----------------------Step %d-----------------------'%step)
        params['step'] = step
        # timer for CP2K
        t1 = time.time()
        if isxTB:
            #CP2K_methods.make_ot_input(params)
            CP2K_methods.run_cp2k_xtb(params)
            molden_filename = F'Diag_{step}-libra-1_0.molden'
        else:
            CP2K_methods.CP2K_input_static( params['cp2k_diag_input_template'], 'Diag_libra', params['trajectory_xyz_filename'], step )
            os.system(F'{mpi_executable} -n {nprocs} {cp2k_exe} -i Diag_libra-{step}.inp -o step_{step}.log')
            molden_filename = F'Diag_libra-{step}-1_0.molden'
        print('Done with step', step,'Elapsed time:',time.time()-t1)

        # now if the counter is equal to zero 
        # just compute the MO overlap of that step.
        if counter == 0:
            t1 = time.time()
            print('Creating shell...')
            shell_1, l_vals = molden_methods.molden_file_to_libint_shell(molden_filename,\
                                                                          is_spherical)
            print('Done with creating shell. Elapsed time:', time.time()-t1)
            t1 = time.time()
            print('Reading energies and eigenvectors....')
            eig_vect_1, energies_1 = molden_methods.eigenvectors_molden(molden_filename,\
                                                                        nbasis(shell_1),l_vals)
            print('Done with reading energies and eigenvectors. Elapsed time:', time.time()-t1)
            print('Computing atomic orbital overlap matrix...')
            t1 = time.time()
            AO_S = compute_overlaps(shell_1,shell_1,nprocs)
            if is_periodic:
                cell = []
                cell.append(params['A_cell_vector'])
                cell.append(params['B_cell_vector'])
                cell.append(params['C_cell_vector'])
                cell = np.array(cell)*units.Angst
                translational_vectors = params['translational_vectors']
                for i1 in range(len(translational_vectors)):
                    translational_vector = np.array(translational_vectors[i1])
                    print(F'Computing the AO overlaps between R({translational_vector[0]},{translational_vector[1]},{translational_vector[2]}) and R(0,0,0)')
                    shell_1p, l_vals = molden_methods.molden_file_to_libint_shell(molden_filename, is_spherical, is_periodic, cell, translational_vector)
                    AO_S += compute_overlaps(shell_1,shell_1p, nprocs)

            print('Done with computing atomic orbital overlaps. Elapsed time:', time.time()-t1)
            t1 = time.time()
            print('Turning the MATRIX to numpy array...')
            AO_S = data_conv.MATRIX2nparray(AO_S)
            #scipy.sparse.save_npz(params['res_dir']+'/AO_S.npz', scipy.sparse.csc_matrix(AO_S))
            print('Done with transforming MATRIX 2 numpy array. Elapsed time:', time.time()-t1)
            lowest_orbital = params['lowest_orbital']
            highest_orbital = params['highest_orbital']
            ## Now, we need to resort the eigenvectors based on the new indices
            print('Resorting eigenvectors elements...')
            t1 = time.time()
            new_indices = CP2K_methods.resort_molog_eigenvectors(l_vals)
            eigenvectors_1 = []

            for j in range(len(eig_vect_1)):
                # the new and sorted eigenvector
                eigenvector_1 = eig_vect_1[j]
                eigenvector_1 = eigenvector_1[new_indices]
                # append it to the eigenvectors list
                eigenvectors_1.append(eigenvector_1)
            eigenvectors_1 = np.array(eigenvectors_1)
            if isUKS:
                ## Alpha spin channel
                #alpha_eigenvectors_1 = eigenvectors_1[0::2]
                #alpha_energies_1 = energies_1[0::2]
                ## Beta spin channel
                #beta_eigenvectors_1 = eigenvectors_1[1::2]
                #beta_energies_1 = energies_1[1::2]
                eig_vec_shape = int(eigenvectors_1.shape[0]/2)
                # Alpha spin channel
                alpha_eigenvectors_1 = eigenvectors_1[0:eig_vec_shape]
                alpha_energies_1 = energies_1[0:eig_vec_shape]
                # Beta spin channel
                beta_eigenvectors_1 = eigenvectors_1[eig_vec_shape:]
                beta_energies_1 = energies_1[eig_vec_shape:]


            print('Done with resorting eigenvectors elements. Elapsed time:',time.time()-t1)
            ##
            t1 = time.time()
            print('Computing and saving molecular orbital overlaps...')
            # Note that we choose the data from lowest_orbital to highest_orbital
            # the values for lowest_orbital and highest_orbital start from 1
            if isUKS:
                S_alpha = np.linalg.multi_dot([alpha_eigenvectors_1, AO_S, alpha_eigenvectors_1.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]
                S_beta  = np.linalg.multi_dot([beta_eigenvectors_1, AO_S, beta_eigenvectors_1.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]
                # creating zero matrix                                     
                mat_block_size = (len(S_alpha), len(S_beta))
                zero_mat = np.zeros(mat_block_size)
                S_step = data_conv.form_block_matrix(S_alpha,zero_mat,zero_mat.T,S_beta)
                # Since a lot of the data are zeros we save them as sparse matrices
                S_step_sparse = scipy.sparse.csc_matrix(S_step)
                E_step = data_conv.form_block_matrix(np.diag(alpha_energies_1)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital],\
                                                     zero_mat,zero_mat.T,\
                                                     np.diag(beta_energies_1)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital])
                E_step_sparse = scipy.sparse.csc_matrix(E_step)

            else:
                S = np.linalg.multi_dot([eigenvectors_1, AO_S, eigenvectors_1.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]
                # creating zero matrix
                mat_block_size = len(S)
                zero_mat = np.zeros((mat_block_size,mat_block_size))
                S_step = data_conv.form_block_matrix(S,zero_mat,zero_mat,S)
                # Since a lot of the data are zeros we save them as sparse matrices
                S_step_sparse = scipy.sparse.csc_matrix(S_step)
                E_step = data_conv.form_block_matrix(np.diag(energies_1)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital],\
                                                     zero_mat,zero_mat,\
                                                     np.diag(energies_1)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital])
                E_step_sparse = scipy.sparse.csc_matrix(E_step)

            scipy.sparse.save_npz(params['res_dir']+F'/S_ks_{step}.npz', S_step_sparse)
            scipy.sparse.save_npz(params['res_dir']+F'/E_ks_{step}.npz', E_step_sparse)
            print('Done with computing molecular orbital overlaps. Elapsed time:', time.time()-t1)
            print(F'Done with step {step}.','Elapsed time:',time.time()-t1_all)
            if params['remove_molden']:
                os.system(F'rm {molden_filename}')

        else:
            # The same procedure as above but now we compute time-overlaps as well
            t1 = time.time()
            print('Creating shell...')
            shell_2, l_vals = molden_methods.molden_file_to_libint_shell(molden_filename,\
                                                                          is_spherical)
            print('Done with creating shell. Elapsed time:', time.time()-t1)
            t1 = time.time()
            print('Reading energies and eigenvectors....')
            eig_vect_2, energies_2 = molden_methods.eigenvectors_molden(molden_filename,\
                                                                        nbasis(shell_2),l_vals)
            print('Done with reading energies and eigenvectors. Elapsed time:', time.time()-t1)
            print('Computing atomic orbital overlap matrix...')
            t1 = time.time()
            AO_S = compute_overlaps(shell_2,shell_2,nprocs)
            AO_St = compute_overlaps(shell_1,shell_2,nprocs)
            if is_periodic:
                cell = []
                cell.append(params['A_cell_vector'])
                cell.append(params['B_cell_vector'])
                cell.append(params['C_cell_vector'])
                cell = np.array(cell)*units.Angst
                translational_vectors = params['translational_vectors']
                for i1 in range(len(translational_vectors)):
                    translational_vector = np.array(translational_vectors[i1])
                    print(F'Computing the AO overlaps between R({translational_vector[0]},{translational_vector[1]},{translational_vector[2]}) and R(0,0,0)')
                    shell_2p, l_vals = molden_methods.molden_file_to_libint_shell(molden_filename, is_spherical, is_periodic, cell, translational_vector)
                    AO_S += compute_overlaps(shell_2,shell_2p, nprocs)
                    AO_St += compute_overlaps(shell_1,shell_2p, nprocs)

            print('Done with computing atomic orbital overlaps. Elapsed time:', time.time()-t1)
            t1 = time.time()
            print('Turning the MATRIX to numpy array...')
            AO_S = data_conv.MATRIX2nparray(AO_S)
            AO_St = data_conv.MATRIX2nparray(AO_St)
            print('Done with transforming MATRIX 2 numpy array. Elapsed time:', time.time()-t1)
            ## Now, we need to resort the eigenvectors based on the new indices
            print('Resorting eigenvectors elements...')
            t1 = time.time()
            new_indices = CP2K_methods.resort_molog_eigenvectors(l_vals)
            eigenvectors_2 = []

            for j in range(len(eig_vect_2)):
                # the new and sorted eigenvector
                eigenvector_2 = eig_vect_2[j]
                eigenvector_2 = eigenvector_2[new_indices]
                # append it to the eigenvectors list
                eigenvectors_2.append(eigenvector_2)
            eigenvectors_2 = np.array(eigenvectors_2)

            if isUKS:
                ## Alpha spin channel
                #alpha_eigenvectors_2 = eigenvectors_2[0::2]
                #alpha_energies_2 = energies_2[0::2]
                ## Beta spin channel
                #beta_eigenvectors_2 = eigenvectors_2[1::2]
                #beta_energies_2 = energies_2[1::2]
                eig_vec_shape = int(eigenvectors_2.shape[0]/2)

                # Alpha spin channel
                alpha_eigenvectors_1 = eigenvectors_1[0:eig_vec_shape]
                alpha_energies_1 = energies_1[0:eig_vec_shape]
                # Beta spin channel
                beta_eigenvectors_1 = eigenvectors_1[eig_vec_shape:]
                beta_energies_1 = energies_1[eig_vec_shape:]

                # Alpha spin channel
                alpha_eigenvectors_2 = eigenvectors_2[0:eig_vec_shape]
                alpha_energies_2 = energies_2[0:eig_vec_shape]
                # Beta spin channel
                beta_eigenvectors_2 = eigenvectors_2[eig_vec_shape:]
                beta_energies_2 = energies_2[eig_vec_shape:]

            print('Done with resorting eigenvectors elements. Elapsed time:',time.time()-t1)
            ##
            t1 = time.time()
            print('Computing and saving molecular orbital overlaps...')

            if isUKS:
                S_alpha = np.linalg.multi_dot([alpha_eigenvectors_2, AO_S, alpha_eigenvectors_2.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]
                St_alpha = np.linalg.multi_dot([alpha_eigenvectors_1, AO_St, alpha_eigenvectors_2.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]

                S_beta = np.linalg.multi_dot([beta_eigenvectors_2, AO_S, beta_eigenvectors_2.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]
                St_beta = np.linalg.multi_dot([beta_eigenvectors_1, AO_St, beta_eigenvectors_2.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]

                S_step = data_conv.form_block_matrix(S_alpha,zero_mat,zero_mat.T,S_beta)
                St_step = data_conv.form_block_matrix(St_alpha,zero_mat,zero_mat.T,St_beta)
                S_step_sparse = scipy.sparse.csc_matrix(S_step)
                St_step_sparse = scipy.sparse.csc_matrix(St_step)

                E_step = data_conv.form_block_matrix(np.diag(alpha_energies_2)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital],\
                                                     zero_mat,zero_mat.T,\
                                                     np.diag(beta_energies_2)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital])

                E_step_sparse = scipy.sparse.csc_matrix(E_step)
            else:
                S = np.linalg.multi_dot([eigenvectors_2, AO_S, eigenvectors_2.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]
                St = np.linalg.multi_dot([eigenvectors_1, AO_St, eigenvectors_2.T])[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital]
                S_step = data_conv.form_block_matrix(S,zero_mat,zero_mat,S)
                St_step = data_conv.form_block_matrix(St,zero_mat,zero_mat,St)
                S_step_sparse = scipy.sparse.csc_matrix(S_step)
                St_step_sparse = scipy.sparse.csc_matrix(St_step)
                E_step = data_conv.form_block_matrix(np.diag(energies_2)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital],\
                                                     zero_mat,zero_mat,\
                                                     np.diag(energies_2)[lowest_orbital-1:highest_orbital,lowest_orbital-1:highest_orbital])
                E_step_sparse = scipy.sparse.csc_matrix(E_step)

            scipy.sparse.save_npz(params['res_dir']+F'/S_ks_{step}.npz', S_step_sparse)
            scipy.sparse.save_npz(params['res_dir']+F'/St_ks_{step-1}.npz', St_step_sparse)
            scipy.sparse.save_npz(params['res_dir']+F'/E_ks_{step}.npz', E_step_sparse)

            print('Done with computing molecular orbital overlaps. Elapsed time:', time.time()-t1)

            shell_1 = shell_2
            energies_1 = energies_2
            eigenvectors_1 = eigenvectors_2

            if params['remove_molden']:
                os.system(F'rm {molden_filename}')

            if isxTB:
                print('Removing unnecessary wfn files...')
                os.system(F'rm OT_{step-1}-RESTART*')
                os.system(F'rm Diag_{step-1}-RESTART*')
            else:
                os.system(F'rm Diag_libra-{step-1}-RESTART*')

            if params['cube_visualization']:
                print('Plotting cube files using VMD...')
                num_orbitals = E_step.shape[0]
                for c_plot, state_to_plot in enumerate(states_to_plot):
                    if isUKS:
                        state_index = state_to_plot-params['lowest_orbital']
                        # Alpha orbitals
                        if isxTB:
                            cube_file_name = F'Diag_{step-1}-WFN_{str(state_to_plot).zfill(5)}_1-1_0.cube'
                        else:
                            cube_file_name = F'Diag_libra-{step-1}-WFN_{str(state_to_plot).zfill(5)}_1-1_0.cube'
                        cube_file_methods.plot_cube_v2(params, cube_file_name, phase_factors_alpha[c_plot,0])
                        # Beta orbitals
                        if isxTB:
                            cube_file_name = F'Diag_{step-1}-WFN_{str(state_to_plot).zfill(5)}_2-1_0.cube'
                        else:
                            cube_file_name = F'Diag_libra-{step-1}-WFN_{str(state_to_plot).zfill(5)}_2-1_0.cube'
                        cube_file_methods.plot_cube_v2(params, cube_file_name, phase_factors_beta[c_plot,0])

                        if params['plot_phase_corrected']:
                            if St_step[state_index, state_index] > 0:
                                f_alpha = 1
                            else:
                                f_alpha = -1

                            if St_step[state_index+num_orbitals, state_index+num_orbitals]>0:
                                f_beta = 1
                            else:
                                f_beta = -1
                            phase_factors_alpha[c_plot,0] *= f_alpha
                            phase_factors_alpha[c_plot,0] *= f_beta

                    else:
                        # Only alpha orbitals
                        state_index = state_to_plot-params['lowest_orbital']
                        if isxTB:
                            cube_file_name = F'Diag_{step-1}-WFN_{str(state_to_plot).zfill(5)}_1-1_0.cube'
                        else:
                            cube_file_name = F'Diag_libra-{step-1}-WFN_{str(state_to_plot).zfill(5)}_1-1_0.cube'
                        cube_file_methods.plot_cube_v2(params, cube_file_name, phase_factors_alpha[c_plot,0])

                        if params['plot_phase_corrected']:
                            if St_step[state_index, state_index] > 0:
                                f_alpha = 1
                            else:
                                f_alpha = -1
                            phase_factors_alpha[c_plot,0] *= f_alpha
                    if params['together_mode']:
                        break
            print(F'Done with step {step}.','Elapsed time:',time.time()-t1_all)
        counter += 1
    # Finally move all the pdos and log files to all_pdosfiles and all_logfiles
    os.system(F'mv *pdos {params["all_pdosfiles"]}/.')
    os.system(F'mv *log {params["all_logfiles"]}/.')
    os.system(F'mv *{params["image_format"].lower()} {params["all_images"]}/.')
    os.system('rm *.wfn* *.tdwfn*')
    if params["remove_cube"]:
        os.system("rm *.cube")
    print('Done with the job!!!')


