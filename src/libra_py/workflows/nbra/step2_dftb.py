# ***********************************************************
# * Copyright (C) 2019-2023 Brendan Smith and Alexey V. Akimov
# * Copyright (C) 2019 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
# ***********************************************************/

"""
.. module:: step2_dftb
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing NAC calculaitons in the 1-electron
       basis of KS orbitals, with the DFTB+ package. It also implements functions for obtaining
       TD-DFTB energies with the DFTB+ package.

.. moduleauthors:: Brendan Smith, Alexey V. Akimov

"""

import os
import sys
import numpy as np

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *

import libra_py.packages.dftbplus.methods as DFTB_methods
from libra_py import units
from libra_py import data_io, eigensolvers, data_conv
import util.libutil as comn


def create_odin_inp(params):
    """
    Creates the input file for the Odin program:
    https://github.com/thomas-niehaus/odin

    Args:
    * params[`filename`] (string) - the name of the .gen file for which we want to compute overlaps 
        [default: "x1.gen"]
    * params[`slakos_prefix`] (string) - the name of the prefix (folder) that contains the Slater-Koster files
        [default: "./"]
    * params[`max_ang_mom`] (dict) - dictionary defining the maximal angular momentum for each element
        [default: {"H":1, "C":2}]
        
    Returns:
    (string): the content of the input file for the ODIN program
    
    """
    filename = params.get("filename", "x1.gen") # gen file
    slakos_prefix = params.get("slakos_prefix", "./") # prefix to the file where the Slater-Koster files are located
    max_l = params.get("max_ang_mom", {"H":1, "C":2} )  # maximal angular momentum numbers for elements

    res = F"'{filename}'\n'{slakos_prefix}'\n'-'\n'.skf'\n"
    f = open(filename)
    A = f.readlines()
    f.close()
    
    elts = A[1].split()
    print(elts)
    for elt in elts:
        if elt in max_l.keys():
            res = F"{res} {max_l[elt]} "
        else:
            print(F"Element: {elt} is not present in the input `max_l` parameter\n max_l = {max_l}")
            sys.exit(0)

    return res

def run_odin(params):
    """
    Run the overlap calculations using ODIN program:
    https://github.com/thomas-niehaus/odin

    Args:
    * params[`ODIN_EXE`] (string) - the path/name of the ODIN executable [default: "odin"]
    * params[`filename`] (string) - the name of the .gen file for which we want to compute overlaps 
        [default: "x1.gen"]
    * params[`slakos_prefix`] (string) - the name of the prefix (folder) that contains the Slater-Koster files
        [default: "./"]
    * params[`max_ang_mom`] (dict) - dictionary defining the maximal angular momentum for each element
        [default: {"H":1, "C":2}]
        
    Returns:
    None: doesn't return anything but executed the ODIN code on the .gen file defined by the `filename`
    to compute the AO overlap matrix `oversqr.dat`
    
    """
    EXE = params.get("ODIN_EXE", "odin")
    inp = create_odin_inp(params)
    f = open("odin.inp", "w")
    f.write(F"{inp}")
    f.close()
    os.system(F"{EXE} < odin.inp")
    


def do_step(snap, params):
    """

    Runs a single-point SCF calculation for a given geometry along a trajectory

    Args:
        snap ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation

            * **params["EXE"]** ( string ): path to the DFTB+ executable [ default: dftb+ ]
            * **params["mo_active_space"]** ( list of ints or None ): indices of the MOs we care about
                The indexing starts from 0, not 1! If set to None - all MOs will be returned. [default: None]
            * **params["md_file"]** ( string ): the name of the xyz file containing the trajectory - the
                file should be in the xyz format produced by the DFTB+ program. [default: "md.xyz"]
            * **params["sp_gen_file"]** ( string ): the name of the .gen file that is listed in the
                DFTB+ input file and contains the geometry of the system (the content of this file will be
                updated for every i value). [default: "x1.gen" ]
            * **params["syst_spec"]** ( string ): the string that is a part of the DFTB+ .gen file and defines
                whether the system is a non-periodic/cluster ("C") or periodic ("S"). [default: "C"]
            * **params["scf_in_file"]** ( string ): the name of the file containing the template for running regular
                SCF calculations for a single-point DFTB+ calculations.

                - It should use the geometry file defined by params["sp_gen_file"].

                [default: "dftb_in_ham1.hsd"]

            * **params["hs_in_file"]** ( string ): the name of the file containing the template for running
                calculations that construct H and S matrices and print them out by the DFTB+ calculations.

                - It should use the geometry file defined by params["sp_gen_file"].
                - It should have the section: "ReadInitialCharges = Yes" to use the previously-converged charge density
                - It should have the section: "WriteHS = Yes" to initialize the writing of the H and S matrices

                [default: "dftb_in_ham2.hsd"]

            * **params["do_tddftb"]** ( bool ): the parameter to control the type of energies to use in the
                calculations. The DFTB+ input files should be setup accordingly

                - True : compute and read TD-DFTB energies. In this case, only the first entry of the returned
                    results ( energies ) is meaningful. Other properties are just None for now

                - False : compute and read only the single-prticle energies [ default ]

            * **params["eigensolver"]** ( string ): the type of the eigensolver to use

               - "libra" : native solver from Libra
               - "eigh" : eigh for Hermitian matrices from scipy [ default ]
               - "eig" : eig for non-Hermitian matrices from scipy
               - "cholesky" : Cholesky-based decomposition from scipy


    Returns:
        tuple: (Ei, MOi, Hi, Si), where:

            * Ei ( CMATRIX(M, M) ), the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us

            * MOi ( CMATRIX(A, M) ), the matrix of the converged Hamiltonian eigenvalues (at given geometry)
                Here, M = len(params["mo_active_space"]) - we output only the MO energies that are of interest to us.
                A - is the number of AOs in this calculation

            * Hi (list of CMATRIX(A, A) ): the Hamiltonian matrices in the AO basis, for each k-point

            * Si (list of CMATRIX(A, A) ): the overlap matrices in the AO basis, for each k-point

    """

    # Now try to get parameters from the input
    critical_params = []
    default_params = {"EXE": "dftb+",
                      "mo_active_space": None,
                      "md_file": "md.xyz", "sp_gen_file": "x1.gen", "syst_spec": "C",
                      "scf_in_file": "dftb_in_ham1.hsd",
                      "hs_in_file": "dftb_in_ham2.hsd",
                      "do_tddftb": False, "use_single_sd_energy": False,
                      "get_dominant_sd_transitions": False,
                      "eigensolver": "eigh" 
                      }

    comn.check_input(params, default_params, critical_params)

    # Get the parameters
    EXE = params["EXE"]
    md_file = params["md_file"]
    sp_gen_file = params["sp_gen_file"]
    syst_spec = params["syst_spec"]
    scf_in_file = params["scf_in_file"]
    hs_in_file = params["hs_in_file"]
    do_tddftb = params["do_tddftb"]
    use_single_sd_energy = params["use_single_sd_energy"]
    get_dominant_sd_transitions = params["get_dominant_sd_transitions"]
    odin_params = params.get("ODIN_PARAMS", {"ODIN_EXE":"odin",
          "filename":"x1.gen", "slakos_prefix":"./",
          "max_ang_mom":{"C":2, "O":2, "H":1, "Ti":3}})
    esolver = params["eigensolver"]

    # Make an input file for SP calculations
    DFTB_methods.xyz_traj2gen_sp(md_file, "x1.gen", snap, syst_spec)

    # Run SCF calculations and generate the charge density for a converged calculations
    # The file x1.gen is used as a geometry
    os.system("cp %s dftb_in.hsd" % scf_in_file)
    os.system("%s" % EXE)

    if do_tddftb:

        # We have chosen to perform a td-dftb calculation for our given system.
        # We need to go into the file EXC.dat file and extract the energies.
        # First, check to see if the EXC.dat file is even there.

        if os.path.isfile('EXC.DAT'):
            print("Found EXC.DAT")
        else:
            print("\nCannot find the file EXC.DAT")
            print(F"Hint: Current working directory is: {os.getcwd()}")
            print("Is this where you expect the file EXC.DAT to be found?")
            print("Check your dftb_in.hsd file to see if excited states calculation is enabled")
            print("Exiting program...\n")
            sys.exit(0)

        configuration_energies = []
        dominant_sd_transitions = []

        f = open("EXC.DAT")
        A = f.readlines()
        f.close()

        sz = len(A)
        for i in range(sz):
            b = A[i].strip().split()

            if not b:
                continue
            else:
                if b[0].replace('.', '', 1).isdigit():

                    if use_single_sd_energy:
                        # We are choosing to use the energy of the dominant SD transition to approximates the
                        # multi-configurational energy. This SD transition is the one with the largest
                        # coefficient in the multi-configurational picture.
                        configuration_energies.append(float(b[6]))

                    else:
                        # Use the multi-configurational energy.
                        configuration_energies.append(float(b[0]))

                    if get_dominant_sd_transitions:
                        # Get the dominant SD transition in either the multi-configurational or single SD picture.
                        dominant_sd_transitions.append([int(b[2]), int(b[4])])

        if use_single_sd_energy:

            # Sort these energies due to the potential for smaller energy values to be higher in index

            # print ("\nIndex and energies of the single SD states BEFORE sorting:\n")
            # print ( [i for i in range(len(configuration_energies))] )
            # print (configuration_energies)
            # print (dominant_sd_transitions)

            single_sd_energies_index_sorted = sorted(
                range(
                    len(configuration_energies)),
                key=lambda k: configuration_energies[k])
            single_sd_energies_sorted = [configuration_energies[i] for i in single_sd_energies_index_sorted]

            if get_dominant_sd_transitions:
                # Sort the dominant SD transitions from those lowest in energy to highest
                dominant_sd_transitions_sorted = [dominant_sd_transitions[i] for i in single_sd_energies_index_sorted]
                with open(params["out_dir"] + "/sd_transitions_%i" % (snap), "w") as out_file:
                    for transition in dominant_sd_transitions_sorted:
                        for index in transition:
                            out_file.write("%s " % (str(index)))

            # print ("\nIndex and energies of the single SD states AFTER sorting:")
            # print (single_sd_energies_index_sorted)
            # print (single_sd_energies_sorted)
            # print (dominant_sd_transitions_sorted)

            num_configs = len(single_sd_energies_sorted)
            E = CMATRIX(num_configs, num_configs)

            for i in range(num_configs):
                E.set(i, i, single_sd_energies_sorted[i])

        else:
            print("Your Multi-Configurational energies are:")
            print(configuration_energies)

            num_configs = len(configuration_energies)
            E = CMATRIX(num_configs, num_configs)

            for i in range(num_configs):
                E.set(i, i, configuration_energies[i])

        # For now, only return E
        return E, None, None, None

    else:

        # Here, we just use the energies of the bands themselves. This corresponds to the single particle picture

        # Just generate the Hamiltonian corresponding to the converged density matrix
        os.system("cp %s dftb_in.hsd" % hs_in_file)
        os.system("%s" % EXE)

        # [0] is because we extract just the gamma-point
        F = DFTB_methods.get_dftb_matrices("hamsqr1.dat")


        #odin_params.update({"filename":"x1.gen"})
        #run_odin(odin_params)
        S = DFTB_methods.get_dftb_matrices("oversqr.dat")


        # Get the dimensions
        ao_sz = F[0].num_of_cols
        ao_act_sp = list(range(0, ao_sz))

        mo_sz = ao_sz
        mo_act_sp = list(range(0, mo_sz))

        if params["mo_active_space"] is not None:
            mo_sz = len(params["mo_active_space"])
            mo_act_sp = list(params["mo_active_space"])

        # Solve the eigenvalue problem with the converged Fock matrix
        # get the converged MOs
        E, MO = None, None
        if esolver == "libra": 
            E = CMATRIX(mo_sz, mo_sz)
            MO = CMATRIX(ao_sz, mo_sz)
            solve_eigen(F[0], S[0], E, MO, 0)

        else:
            np_F = data_conv.MATRIX2nparray(F[0])
            np_S = data_conv.MATRIX2nparray(S[0])
        
            if esolver == "eigh":
                eigvals, eigvecs = eigensolvers.generalized_eigensolve_scipy(0.5*(np_F + np_F.T), 0.5*(np_S+np_S.T), hermitian=True, sort=True)
            elif esolver == "eig":
                eigvals, eigvecs = eigensolvers.generalized_eigensolve_scipy(np_F, np_S, hermitian=False, sort=True)
            elif esolver == "cholesky":
                eigvals, eigvecs = eigensolvers.generalized_eigensolve_scipy_cholesky(0.5*(np_F + np_F.T), 0.5*(np_S+np_S.T), hermitian=True)

            E = data_conv.nparray2CMATRIX(np.diag(eigvals) )
            MO = data_conv.nparray2CMATRIX(eigvecs)

        print(F" is Fock matrix symmetric?  { (F[0] - F[0].H()).real().max_elt() }")
        print(F" is ovlp matrix symmetric?  { (S[0] - S[0].H()).real().max_elt() }") 
        print(F"Solving eigenvalue problem - max elt of error =  {(F[0] * MO - S[0] * MO * E).real().max_elt()}")

        # Extract the E sub-matrix
        E_sub = CMATRIX(mo_sz, mo_sz)
        pop_submatrix(E, E_sub, mo_act_sp, mo_act_sp)

        # Extract the MO sub-matrix
        MO_sub = CMATRIX(ao_sz, mo_sz)
        pop_submatrix(MO, MO_sub, ao_act_sp, mo_act_sp)

        print(F"Checking the sub-eigenvalue problem - max elt of error =  {(F[0] * MO_sub - S[0] * MO_sub * E_sub).real().max_elt()}")
        iden = MATRIX(ao_sz, ao_sz); iden.identity();
        print(F"Checking the orthogonality { ((MO.H() * S[0] * MO).real() - iden).max_elt() }")

        return E_sub, MO_sub, F, S


def do_ovlp(snap, params):
    """

    Compute the overlap matrix in the AO basis for two geometries, i and i+1

    Args:
        snap ( int ): index of the time step to be used from the trajectory file
        params ( dictionary ): the control parameters of the simulation

            * **params["EXE"]** ( string ): path to the DFTB+ executable [ default: dftb+ ]
            * **params["ODIN_PARAMS"]** ( dict ): dictionary containing parameters of the `run_odin` function
            * **params["md_file"]** ( string ): the name of the xyz file containing the trajectory - the
                file should be in the xyz format produced by the DFTB+ program. [default: "md.xyz"]
            * **params["ovlp_gen_file"]** ( string ): the name of the .gen file that is listed in the
                DFTB+ input file and contains the geometry of the system at two time steps
                (the content of this file will be updated for every i value). [default: "x2.gen" ]
            * **params["syst_spec"]** ( string ): the string that is a part of the DFTB+ .gen file and defines
                whether the system is a non-periodic/cluster ("C") or periodic ("S"). [default: "C"]
            * **params["ovlp_in_file"]** ( string ): the name of the file containing the template for running
                calculations that construct H and S matrices and print them out by the DFTB+ calculations.

                - It should use the geometry file defined by params["ovlp_gen_file"].
                - It should have the section: "WriteHS = Yes" to initialize the writing of the H and S matrices

                [default: "dftb_in_overlaps.hsd"]

            * **params["tol"]** (float ): the threshold for fixing the elements of the off-diagonal blocks. The off-diagonal block elements
                are only reset to the corresponding values of the on-diagonal blocks (or to zeros) if they exceed them by this
                amount. The larger values of this parameter will favor keeping the off-diagonal block matrix elements to
                what they are in the direct calculations, which may be problematic for very close distances (e.g. identical geometries).
                On the contrary, the smaller values of this parameter will favor replacing the off-diagonal block matrix elements
                with the corresponding on-diagonal block matrix elements or zeroes. This may artificially decrease the NACs and may
                slow down the nonadiabatic transitions. E.g. in the extreme case of `tol = 0.0`, all NACs would be zero, so there will
                be no dynamics. [ default: 0.5 ]


    Returns:
        CMATRIX(A, A): the matrix of the AO overlaps for two geometries, where A - is the size of the AO basis

    """

    # Now try to get parameters from the input
    critical_params = []
    default_params = {"EXE": "dftb+", "ODIN_EXE":"odin",
                      "md_file": "md.xyz", "ovlp_gen_file": "x2.gen", "syst_spec": "C",
                      "ovlp_in_file": "dftb_in_overlaps.hsd", "tol": 2.5
                      }
    comn.check_input(params, default_params, critical_params)

    # Get the parameters
    EXE = params["EXE"]
    md_file = params["md_file"]
    ovlp_gen_file = params["ovlp_gen_file"]
    syst_spec = params["syst_spec"]
    ovlp_in_file = params["ovlp_in_file"]
    tol = params["tol"]

    odin_params = params.get("ODIN_PARAMS", {"ODIN_EXE":"odin",
          "filename":"x2.gen", "slakos_prefix":"./",
          "max_ang_mom":{"C":2, "O":2, "H":1, "Ti":3}})

    # Make an input file for the overlap calculations
    DFTB_methods.xyz_traj2gen_ovlp(md_file, "x2.gen", snap, snap + 1, syst_spec)

    # Run SCF calculations and generate the charge density for a converged calculations
    # The file x2.gen is used as a geometry
    #os.system("cp %s dftb_in.hsd" % ovlp_in_file)
    #os.system("%s" % EXE)
    odin_params.update({"filename":"x2.gen"})
    run_odin(odin_params)

    # Get the Hamiltonian
    Sbig = DFTB_methods.get_dftb_matrices("oversqr.dat") #, None, None, 0, tol)
    norbs = int(Sbig[0].num_of_cols / 2)

    act_sp1 = list(range(0, norbs))
    act_sp2 = list(range(norbs, 2 * norbs))
    S = CMATRIX(norbs, norbs)
    pop_submatrix(Sbig[0], S, act_sp1, act_sp2)

    S00 = CMATRIX(norbs, norbs)
    pop_submatrix(Sbig[0], S00, act_sp1, act_sp1)

    S11 = CMATRIX(norbs, norbs)
    pop_submatrix(Sbig[0], S11, act_sp2, act_sp2)

    return S, S00, S11


def run_step2(params):
    """

    Calculate the overlaps, transition dipole moments, and vibronic Hamiltonian matrix elements in the AO basis

    Currently under development 1/30/2020 - This function is currently not compatable with the parameters
                                            "get_midpoint_energy" or "do_tddftb".

    Args:
        params ( dictionary ): the control parameters of the simulation

            * **params["dt"]** ( double ): nuclear dynamics timestep - as encoded in the trajectory [ units: a.u., default: 41.0 ]
            * **params["isnap"]** ( int ): initial frame  [ default: 0 ]
            * **params["fsnap"]** ( int ): final frame  [ default: 1 ]
            * **params["out_dir"]** ( string ): the path to the directory that will collect all the results
                If the directory doesn't exist, it will be created  [ default: "res" ]

            SeeAlso:  the description of the parameters in ```do_ovlp(i, params)``` and in ```do_step(i, params)```

    Returns:
        None: but generates the files (S, St, and hvib) indexed in the range [isnap, fsnap), for
        instance, if isnap = 0, fsnap = 3, we will have files hvib_0, hvib_1, hvib_2


    """

    critical_params = []
    default_params = {"dt": 1.0 * units.fs2au, "isnap": 0, "fsnap": 1, "out_dir": "res"}
    comn.check_input(params, default_params, critical_params)

    # Get the parameters
    dt = params["dt"]
    isnap = params["isnap"]
    fsnap = params["fsnap"]
    out_dir = params["out_dir"]

    # Create <out_dir> directory if it does not exist yet
    if os.path.isdir(out_dir):
        pass
    else:
        os.system("mkdir %s" % out_dir)

    print("In run_step2 : params = ", params )

    # Compute
    E_prev, U_prev, Hao_prev, Sao_prev = do_step(isnap, params)
    #E_prev, U_prev, Hao_prev, Sao_prev = do_step(0, params) # just a test

    os.system(F"cp hamsqr1.dat {out_dir}/hamsqr1_step{isnap}.dat")
    os.system(F"cp oversqr.dat {out_dir}/oversqr_small_step{isnap}.dat")

    for i in range(isnap + 1, fsnap - 1):
        E_curr, U_curr, Hao_curr, Sao_curr = do_step(i, params)
        #E_curr, U_curr, Hao_curr, Sao_curr = do_step(0, params)  # jut a test
        #print("MO dimensions = ", U_curr.num_of_rows, U_curr.num_of_cols)
        #U_curr.show_matrix()

        os.system(F"cp hamsqr1.dat {out_dir}/hamsqr1_step{i}.dat")
        os.system(F"cp oversqr.dat {out_dir}/oversqr_small_step{i}.dat")

        S, S00, S11 = do_ovlp(i - 1, params)
        #S, S00, S11 = do_ovlp(0, params)  # just a test

        os.system(F"cp oversqr.dat {out_dir}/oversqr_big_step{isnap}.dat")

        print(F"S00 - Sao_prev = {(S00 - Sao_prev[0]).real().max_elt() }") 

        # S.real().show_matrix("res/AOS_%i_re" % (i) )

        TDM = U_prev.H() * S * U_curr
        #TDM = U_prev.H() * S00 * U_curr  # just a test
        #print("TDM", TDM.num_of_rows, TDM.num_of_cols)
        #TDM.show_matrix()
        Hvib = 0.5 * (E_prev + E_curr) - (0.5j / dt) * (TDM - TDM.H())

        # Overlaps
        s = 0.5 * (U_prev.H() * S00 * U_prev + U_curr.H() * S11 * U_curr)
        s.real().show_matrix("%s/S_%i_re" % (out_dir, i - 1))
        s.imag().show_matrix("%s/S_%i_im" % (out_dir, i - 1))

        # Time-overlaps
        TDM.real().show_matrix("%s/St_%i_re" % (out_dir, i - 1))
        TDM.imag().show_matrix("%s/St_%i_im" % (out_dir, i - 1))

        # Vibronic Hamiltonians
        Hvib.real().show_matrix("%s/hvib_%i_re" % (out_dir, i - 1))
        Hvib.imag().show_matrix("%s/hvib_%i_im" % (out_dir, i - 1))

        # Current becomes the old
        E_prev = CMATRIX(E_curr)
        U_prev = CMATRIX(U_curr)
        Sao_prev[0] = CMATRIX(Sao_curr[0])


def run_step2_lz(params):
    """

    Calculate the diagonal vibronic Hamiltonian matrix elements (energies) in the AO basis and (optionally) the TD-DFTB energies

    Args:
        params ( dictionary ): the control parameters of the simulation

            * **params["dt"]** ( double ): nuclear dynamics timestep - as encoded in the trajectory [ units: a.u., default: 41.0 ]
            * **params["isnap"]** ( int ): initial frame  [ default: 0 ]
            * **params["fsnap"]** ( int ): final frame  [ default: 1 ]
            * **params["out_dir"]** ( string ): the path to the directory that will collect all the results
                If the directory doesn't exist, it will be created  [ default: "res" ]
            * **params["get_midpoint_energy"]** ( bool ): if True, compute and read energies as Hvib = 0.5*(E_prev + E_curr) (As is done in Pyxaid NBRA)
                                                   if False, compute and read energies at everytime timestep  Hvib = E_curr

            SeeAlso:  the description of the parameters in ```do_step(i, params)```

    Returns:
        None: but generates the files (hvib) indexed in the range [isnap, fsnap), for
        instance, if isnap = 0, fsnap = 3, we will have files hvib_0, hvib_1, hvib_2.
        Only diagonal elements of the matrix are populated.

    """

    critical_params = []
    default_params = {"dt": 1.0 * units.fs2au, "isnap": 0, "fsnap": 1, "out_dir": "res", "get_midpoint_energy": True}
    comn.check_input(params, default_params, critical_params)

    # Get the parameters
    dt = params["dt"]
    isnap = params["isnap"]
    fsnap = params["fsnap"]
    out_dir = params["out_dir"]
    get_midpoint_energy = params["get_midpoint_energy"]

    # Create <out_dir> directory if it does not exist yet
    if os.path.isdir(out_dir):
        pass
    else:
        os.system("mkdir %s" % out_dir)

    if get_midpoint_energy:

        E_prev, U_prev, Hao_prev, Sao_prev = do_step(isnap, params)
        for i in range(isnap + 1, fsnap):

            E_curr, U_curr, Hao_curr, Sao_curr = do_step(i, params)
            Hvib = 0.5 * (E_prev + E_curr)

            Hvib = E_curr
            # Vibronic Hamiltonians - Diag. elements only
            Hvib.real().show_matrix("%s/hvib_%i_re" % (out_dir, i - 1))

            # Current becomes the old
            E_prev = CMATRIX(E_curr)

    else:
        for i in range(isnap, fsnap - 1):

            E_curr, U_curr, Hao_curr, Sao_curr = do_step(i, params)
            Hvib = E_curr
            # Vibronic Hamiltonians - Diag. elements only
            Hvib.real().show_matrix("%s/hvib_%i_re" % (out_dir, i))
