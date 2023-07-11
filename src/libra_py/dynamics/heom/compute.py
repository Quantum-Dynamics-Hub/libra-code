#*********************************************************************************
#* Copyright (C) 2019-2020 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#***********************************************************************************
"""
.. module:: compute
   :platform: Unix, Windows
   :synopsis: This module implements a wrapper function for doing HEOM dynamics
       The code is a translation/refactoring of the Fortran code or Amber Jain & Joe Subotnik
       https://github.com/subotnikgroup/HEOM_Amber

       List of functions:
           * run_dynamics(dyn_params, Ham, rho_init)

.. moduleauthors:: Alexey V. Akimov



"""

__author__ = "Alexey V. Akimov"
__copyright__ = "Copyright 2019-200 Alexey V. Akimov"
__credits__ = ["Amber Jain, Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"


import os
import sys
import math
import copy
import time

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


import util.libutil as comn
import libra_py.units as units
from . import save

def aux_print_matrices(step, x):
    print(F"= step = {step} =")
    nmat = len(x)
    nrows = x[0].num_of_rows
    ncols = x[0].num_of_cols
    for imat in range(nmat):
        print(F"== imat = {imat} ==")
        for row in range(nrows):
            line = ""
            for col in range(ncols):
                line = line + F"{x[imat].get(row, col)}  "
            print(line)


def update_filters(rho_scaled, params, aux_memory):
    """

    This function takes the current ADMs hierarchy, computes the corresponding time-derivatives
    (based on the input params["adm_list"]), then it updates the params["adm_list"] to indicate
    the equations for which the |dADM/dt (max element)| is larger than a given threshold,
    params["adm_deriv_tolerance"].

    This function also runs over all ADMs to determine those which have |ADM (max element)| larger
    than the params["adm_tolerance"]. This will determine the update of the params["nonzero"] lists.
    The ADMs that are determined to be "zero" (discarded until later) can also be set to 0.0, the
    option controlled by params["do_zeroing"].
    In the end, the variable `rho` is updated with the correspondingly updated variables

    aux_memory should have allocated:
    - rho_unpacked_scaled
    - drho_unpacked_scaled

    """


    unpack_mtx(aux_memory["rho_unpacked_scaled"], rho_scaled)

    drho_scaled = compute_heom_derivatives(rho_scaled, params)
    unpack_mtx(aux_memory["drho_unpacked_scaled"], drho_scaled)

    # Filtering of the derivatives - defines the active EOMs, params["adm_list"]
    trash = filter(aux_memory["drho_unpacked_scaled"], params["adm_list"], params["adm_deriv_tolerance"], 0)

    # Filtering of the densities - defines the list of nonzero derivatives, params["nonzero"]
    params["nonzero"] = filter(aux_memory["rho_unpacked_scaled"], trash, params["adm_tolerance"], params["do_zeroing"])

    pack_mtx(aux_memory["rho_unpacked_scaled"], rho_scaled)





def transform_adm(rho, rho_scaled, aux_memory, params, direction):
    """
    aux_memory should have allocated:
    - rho_unpacked
    - rho_unpacked_scaled

    direction:
     1 :  raw -> scaled
    -1 :  scaled -> raw

    """

    nn_tot = len(aux_memory["rho_unpacked"])

    # Forward
    if direction == 1:

        unpack_mtx(aux_memory["rho_unpacked"], rho)

        if params["do_scale"]==0:
            for n in range(nn_tot):
                aux_memory["rho_unpacked_scaled"][n] = CMATRIX(aux_memory["rho_unpacked"][n])

        elif params["do_scale"]==1:
            scale_rho(aux_memory["rho_unpacked"], aux_memory["rho_unpacked_scaled"], params)

        pack_mtx(aux_memory["rho_unpacked_scaled"], rho_scaled)


    # Backward
    elif direction == -1:

        unpack_mtx(aux_memory["rho_unpacked_scaled"], rho_scaled)

        # We want to save only actual ADMs so we need to convert the
        # scaled one back to the unscaled
        if params["do_scale"]==0:
            for n in range(nn_tot):
                aux_memory["rho_unpacked"][n] = CMATRIX(aux_memory["rho_unpacked_scaled"][n])

        elif params["do_scale"]==1:
            # rho_unpacked_scaled -> rho_unpacked
            inv_scale_rho(aux_memory["rho_unpacked"], aux_memory["rho_unpacked_scaled"], params)

        pack_mtx(aux_memory["rho_unpacked"], rho)


def run_dynamics(dyn_params, Ham, rho_init):

    """
    This functions integrates the HEOM for a given system's Hamiltonian, initial conditions, and bath parameters

    Args:
        dyn_params ( dictionary )
            Parameters controlling the execution of the dynamics. Can contain:

            =============== Properties of the bath ================

            * **dyn_params["KK"]** ( int )
                Defines the the number of Matsubara modes (KK+1) - one needs to
                achieve a convergence w.r.t. this parameter [ default: 0]

            * **dyn_params["LL"]** ( int )
                The hierarchy level - one needs to achieve a convergence w.r.t. this parameter [ default: 10 ]

            * **dyn_params["gamma"]** ( float )
                The system-bath interaction ("collision") frequency, related
                to the bath's "friction" on the system [ units: Ha,  default: 1.0/(0.1 * units.ps2au) ]

            * **dyn_params["eta"]** ( float )
                Reorganization energy of bath [ units: Ha, default: 2.0 * 50.0 * units.inv_cm2Ha ]

            * **dyn_params["temperature"]** ( float )
                Temperature of the bath [ units: K, default: 300 ]

            * **dyn_params["el_phon_couplings"]** ( list of CMATRIX; (nstates+1) x CMATRIX(nstates, nstates) )
                The matrices that describe how each electronic phonon is coupled to various electronic states
                in the simplest picture, when a phonon k is coupled to electronic state m, the matrix el_phon_couplings[k]
                will contain 1.0 at the element (m,m) and zeroes everywhere else. You can of course define more general
                situations when one phonon is coupled to many states (and perhaps to their coherences).
                The convention is: dyn_params["el_phon_couplings"][0] - contains trash, and only the next element is the
                first actual coupling. [ default: see the simplest example described above ]


            =============== Parameters of the dynamics ================

            * **dyn_params["dt"]** ( float )
                Time-integration timestep [ units: a.u. of time, default: 0.1*units.fs2au ]

            * **dyn_params["nsteps"]** ( int )
                How many steps of the dynamics to perform [ default: 10 ]

            * **dyn_params["verbosity"]** ( int )
                The level of the run-time printout of any useful information [ default: -1 ]

            * **dyn_params["progress_frequency"]** ( float )
                The frequency (defined as the as fraction of the entire simulation run)
                to print out the progress message. This number shall be in the interval [0, 1]
                [ default: 0.1 ]


            =============== Algorithmic parameters ================

            * **dyn_params["truncation_scheme"]** ( int )
                How to truncate the HEOM equations. Options are:

                - 0 : no truncation
                - 1 : according to Schulten with real part of Matsubara terms [ default ]
                - 2 : according to Schulten with full Matsubara term
                - 3 : according to Shi with real part of Matsubara terms
                - 4 : according to Shi with full Matsubara terms

            * **dyn_params["do_scale"]** ( int )
                Whether to use the scaled HEOM version. Options are:

                - 0 : don't use the scaled HEOM
                - 1 : use it according to JCP 130, 084105, 2009 [ default ]

            * **dyn_params["tolerance"]** ( float )
                The threshold for discarding the auxiliary density matrices.
                The larger it is, the fewer auxiliary density matrices survives in the HEOM
                and the faster the calculations become. [ default: 1e-6 ]

            * **dyn_params["filter_after_steps"]** ( int )
                Denfines the frequency (after how may steps of the dynamics)
                the auxiliary density matrices will be checked to be potentially discarded
                (based on the `tolerance` parameter set above) [ default: 10 ]

            * **dyn_params["num_threads"]** ( int )
                The number of OMP threads to use to parallelize the calculations [ default: 1 ]

            =============== Parameters for saving data ================

            * **dyn_params["prefix"]** ( string )
                The name of the folder which will contain the results computed.
                This folder will be created if not existent already [ default: "out" ]

            * **dyn_params["hdf5_output_level"]** ( int )
                The level of on-the-fly HDF5 file printout.
                The file is called "data.hdf" and is stored in the directory defined by the `prefix` variable.
                The larger it is, the more properties will be saved, and the larger the size of the generated files will be.
                In this cases, the files are being written as the calculations go, so it may be a quite slow process (the main
                bottleneck of the calculations!). However, if your calculations are stopped before they reach the end,
                you shall still have the data. [ default: 0 ]

            * **dyn_params["txt_output_level"]** ( int )
                The level of the text-based output of the computed results. Is not yet implemented. [ default: 0 ]

            * **dyn_params["mem_output_level"]** ( int )
                The level of the memory-based HDF5 file creation. This is a flag similar to the `hdf5_output_level`
                The file is called "mem_data.hdf" and is stored in the directory defined by the `prefix` variable.
                The larger it is, the more properties will be saved, and the larger the size of the generated files will be.
                Unlike with the "hdf5_output_level", the files are being written only when the calculations end. This way,
                the writing of the files is much much faster than with "hdf5_output_lelel", but if the calculations are stopped
                or crash the mid-way, you will have nothing stored. Also, in this case, the function also returns a tuple with
                all the variables stored, so ready to use in the Python program. Keep in mind that although this option
                is much faster than "hdf5_output_level", all the results are first stored in the OS memory before they are
                dumped into the HDF5 files, so for large systems/calculations you may need good amount of RAM. [ default: 3 ]

            * **dyn_params["properties_to_save"]** ( list of strings )
                The names of the datasets (data) to store to the HDF5 files.
                Note that one needs to satisfy both the *_output_level and list the dataset in this parameter
                in order to have it actually saved in the file.

                Available names:

                - **timestep** ( int )
                    The index of the timestep
                    *_output_level >= 1

                - **time** ( float )
                    The time of dynamics [ in a.u. of time]
                    *_output_level >= 1

                - **denmat** ( CMATRIX(nstates, nstates) )
                    The density matrix evolution
                    *_output_level >= 3

                [ default: [ "timestep", "time", "denmat" ] ]


            * **dyn_params["use_compression"]** ( int )
                Whether to use the data compression (via gzip) when storing data to HDF5 files.
                Options:

                - 0 : don't compress data [ default ]
                - 1 : compress data

                The experience shows that it is better not to use the compression. It gets slower.

            * **dyn_params["compression_level"]** ( 3 x int )
                The level of compression for integers, float, and complex numbers.

                Must be a number between 0 and 9, including the ends.

                [ default: [0,0,0] ]



        Ham ( CMATRIX(nstates, nstates) )
            Define the system's electronic Hamiltonian. The diagonal elements contain
            the site energies, the off-diagonal elements contain electronic couplings [ units: Ha ]

        rho ( CMATRIX(nstates, nstates) )
            Is the initial density matrix describing the quantum system.
            It's dimensions must be the same as those of the `Ham` variable.



    """


    params = dict(dyn_params)


    # Parameters and dimensions
    critical_params = [  ]
    default_params = { "KK":0, "LL":10,
                       "gamma": 1.0/(0.1 * units.ps2au),
                       "eta": 2.0 * 50.0 * units.inv_cm2Ha,
                       "temperature": 300.0,
                       "el_phon_couplings":initialize_el_phonon_couplings(Ham.num_of_cols),

                       "dt":0.1*units.fs2au, "nsteps":10,
                       "verbosity":-1, "progress_frequency":0.1,

                       "truncation_scheme":1, "do_scale":0,
                       "adm_tolerance":1e-6,  "adm_deriv_tolerance":1e-12,
                       "filter_after_steps":1, "do_zeroing":0,
                       "num_threads":1,

                       "prefix":"out",
                       "hdf5_output_level":0, "txt_output_level":0, "mem_output_level":3,
                       "properties_to_save": [ "timestep", "time", "denmat"],
                       "use_compression":0, "compression_level":[0,0,0]
                     }

    comn.check_input(params, default_params, critical_params)

    nsteps = params["nsteps"]
    print_freq = int(params["progress_frequency"]*nsteps)


    #============= System ======================
    params.update({"Ham" : Ham})
    nquant = Ham.num_of_cols


    #============== HEOM topology ==============


    KK = dyn_params["KK"]
    LL = dyn_params["LL"]

    all_vectors = intList2()
    vec_plus = intList2()
    vec_minus = intList2()

    gen_hierarchy(nquant * (KK+1), LL, params["verbosity"], all_vectors, vec_plus, vec_minus)
    params.update( { "nvec":all_vectors, "nvec_plus":vec_plus, "nvec_minus":vec_minus } )

    nn_tot = len(all_vectors)


    all_indices = []
    init_nonzero = []
    for n in range(nn_tot):
        all_indices.append(n)
        init_nonzero.append(1)


    #============ Bath update =====================
    gamma_matsubara = doubleList()
    c_matsubara = complexList()

    setup_bath(KK, params["eta"], params["gamma"], params["temperature"], gamma_matsubara, c_matsubara)
    params.update({ "gamma_matsubara": gamma_matsubara, "c_matsubara":c_matsubara  } )

    if params["verbosity"]>=1:
        for k in range(KK+1):
            print(F" k = {k} gamma_matsubara[{k}] = {gamma_matsubara[k]}  c_matsubara[{k}] = {c_matsubara[k]}")

    #============= Initialization ============

    rho = CMATRIX((nn_tot)*nquant, nquant)  # all rho matrices stacked on top of each other
    rho_scaled = CMATRIX((nn_tot)*nquant, nquant)  # all rho matrices stacked on top of each other
    #drho = CMATRIX((nn_tot)*nquant, nquant)  # all rho matrices stacked on top of each other

    aux_memory = {"rho_unpacked" : CMATRIXList(),
                  "rho_unpacked_scaled" : CMATRIXList(),
                  "drho_unpacked" : CMATRIXList(),
                  "drho_unpacked_scaled" : CMATRIXList()
                 }
    for n in range(nn_tot):
        aux_memory["rho_unpacked"].append( CMATRIX(nquant, nquant))
        aux_memory["rho_unpacked_scaled"].append( CMATRIX(nquant, nquant))
        aux_memory["drho_unpacked"].append( CMATRIX(nquant, nquant))
        aux_memory["drho_unpacked_scaled"].append( CMATRIX(nquant, nquant))

    # Initial conditions
    x_ = Py2Cpp_int(list(range(nquant)))
    y_ = Py2Cpp_int(list(range(nquant)))
    push_submatrix(rho, rho_init, x_, y_)

    #unpack_mtx(aux_memory["rho_unpacked"], rho)


    #========== Scale working ADMs ====================
    if params["verbosity"]>=2 and params["do_scale"]==1:
        print("Scaling factors")
        for n in range(nn_tot):
            for m in range(nquant):
                for k in range(KK+1):
                    n_mk = all_vectors[n][m*(KK+1)+k]
                    scl = 1.0/math.sqrt( FACTORIAL(n_mk) * FAST_POW(abs(c_matsubara[k]), n_mk))
                    print(F" n={n} m={m} k={k}  scaling_factor={scl}")

    # raw -> scaled
    transform_adm(rho, rho_scaled, aux_memory, params, 1)

    #========== Filter scaled ADMs ======================
    params.update({ "nonzero": Py2Cpp_int(init_nonzero), "adm_list": Py2Cpp_int( all_indices ) } )
    update_filters(rho_scaled, params, aux_memory)



    if params["verbosity"]>=2:
        print("nonzero = ", Cpp2Py(params["nonzero"]))
        print("adm_list = ", Cpp2Py(params["adm_list"]))

        if params["verbosity"]>=4:
            print("ADMs")
            aux_print_matrices(0, aux_memory["rho_unpacked"])
            print("Scaled ADMs")
            aux_print_matrices(0, aux_memory["rho_unpacked_scaled"])


    # Initialize savers
    _savers = save.init_heom_savers(params, nquant)

    #============== Propagation =============


    start = time.time()
    for step in range(params["nsteps"]):

        #================ Saving and printout ===================
        # scaled -> raw
        transform_adm(rho, rho_scaled, aux_memory, params, -1)

        # Save the variables
        save.save_heom_data(_savers, step, print_freq, params, aux_memory["rho_unpacked"])


        if step%print_freq==0:
            print(F" step= {step}")

            if params["verbosity"]>=3:
                print("nonzero = ", Cpp2Py(params["nonzero"]))
                print("adm_list = ", Cpp2Py(params["adm_list"]))

            if params["verbosity"]>=4:
                print("ADMs")
                aux_print_matrices(0, aux_memory["rho_unpacked"])
                print("Scaled ADMs")
                aux_print_matrices(0, aux_memory["rho_unpacked_scaled"])




        #============== Update the list of active equations = Filtering  ============
        if step % params["filter_after_steps"] == 0:

            # To assess which equations to discard, lets estimate the time-derivatives of rho
            # for all the matrices
            params["adm_list"] = Py2Cpp_int( all_indices )

            update_filters(rho_scaled, params, aux_memory)


        #================= Propagation for one timestep ==================================
        rho_scaled = RK4(rho_scaled, params["dt"], compute_heom_derivatives, params)


    end = time.time()
    print(F"Calculations took {end - start} seconds")


    # For the mem_saver - store all the results into HDF5 format only at the end of the simulation
    if _savers["mem_saver"] != None:
        prefix = params["prefix"]
        _savers["mem_saver"].save_data( F"{prefix}/mem_data.hdf", params["properties_to_save"], "w")
        return _savers["mem_saver"]
