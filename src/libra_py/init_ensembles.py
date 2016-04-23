#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file init_ensembles.py
# This module implements the functions that allocate memory for various objects

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def init_ext_hamiltonians(ntraj, nnucl, nel, verbose=0):
    ##
    # This function allocates memory for matrices (external objects). The function
    # also creates an array of External_Hamiltonian objects and bind them to the external matrices
    # 
    # \param[in] ntraj The number of trajectories in the ensemble
    # \param[in] nnucl The number of nuclear DOF ( usually, = 3* Number_of_atoms)
    # \param[in] nel The number of electronic DOF
    # \param[in] verbose The parameter controlling the amount of debug printing: 0 (default)
    # - no printing, 1 - print info

    # Create an list of External Hamiltonian objects
    ham = []
    for i in xrange(ntraj):
        ham_i = Hamiltonian_Extern(nel, nnucl) 
        ham_i.set_rep(1)            # adiabatic representation
        ham_i.set_adiabatic_opt(0)  # use the externally-computed 
                                    # adiabatic electronic Hamiltonian and derivatives
        ham_i.set_vibronic_opt(0)   # use the externally-computed vibronic Hamiltonian and derivatives
        ham.append(ham_i)

    # bind actual matrices to the External_Hamiltonian objects
    ham_adi = []; d1ham_adi = []; ham_vib = []


    for i in xrange(ntraj):
        # create external matrix for adiabatic Hamiltonian
        ham_adi_i = MATRIX(nel, nel)
        ham_adi.append(ham_adi_i)
        # bind adiabatic Hamiltonian
        ham[i].bind_ham_adi(ham_adi[i])

        # create external matrices for the derivatives
        d1ham_adi_i= MATRIXList()
        for k in xrange(nnucl):
            tmp = MATRIX(nel, nel)
            d1ham_adi_i.append(tmp)
        d1ham_adi.append(d1ham_adi_i)
        # bind derivative of adiabatic Hamiltonian
        ham[i].bind_d1ham_adi(d1ham_adi[i])

        # create external matrix for vibronic Hamiltonian
        ham_vib_i = CMATRIX(nel, nel)
        ham_vib.append(ham_vib_i)
        # bind vibronic Hamiltonian
        ham[i].bind_ham_vib(ham_vib[i])

    if verbose==1:
        print "Addresses of the working matrices"
        for i in xrange(ntraj):
            print "%i th Hamiltonian is" % i
            print "ham_adi = ", ham_adi[i]
            print "d1ham_adi = ", d1ham_adi[i]
            for k in xrange(nnucl):
                print "d1ham_adi[",k,"]= ", d1ham_adi[i][k]
            print "ham_vib = ", ham_vib[i]

    # Returned values:
    # ham - list of External_Hamiltonian objects
    # ham_adi - list of corresponding matrices
    # d1ham_adi - list of the lists containing derivative matrices
    # ham_vib - list of matrices containing vibronic Hamiltonians 
    return ham, ham_adi, d1ham_adi, ham_vib



def init_mols(syst, ntraj, nnucl, verbose=0):
    ##
    # This function creates a list of Nuclear objects and initializes them by extracting info from the
    # provided systems - in "syst"
    # 
    # \param[in] syst The list of chemsystem objects with the information about nuclear DOF
    # \param[in] ntraj The number of trajectories in the ensemble
    # \param[in] nnucl The number of nuclear DOF ( usually, = 3* Number_of_atoms)
    # \param[in] verbose The parameter controlling the amount of debug printing: 0 (default)
    # - no printing, 1 - print info

    # Initialize nuclear variables
    mol = []
    for i in xrange(ntraj):
        mol.append(Nuclear(nnucl))
        syst[i].extract_atomic_q(mol[i].q)
        syst[i].extract_atomic_p(mol[i].p)
        syst[i].extract_atomic_f(mol[i].f)
        syst[i].extract_atomic_mass(mol[i].mass)

        if verbose==1:
            print i,"mol=",mol[i]
            for k in xrange(syst[i].Number_of_atoms):
                print "mol[%d][%d]=%f,%f,%f"%(i,k,mol[i].q[3*k+0],mol[i].q[3*k+1],mol[i].q[3*k+2])

    # Returned values:
    # mol - the list of the "Nuclear" objects for given set of systems
    return mol


def init_therms(ntraj, nnucl, params, verbose=0):
    ##
    # This function creates a list of Thermostat objects and initializes according to the given parameters
    #
    # \param[in] ntraj The number of trajectories in the ensemble
    # \param[in] nnucl The number of nuclear DOF ( usually, = 3* Number_of_atoms)
    # \param[in] verbose The parameter controlling the amount of debug printing: 0 (default)
    # - no printing, 1 - print info


    # Initialize Thermostat object    
    therm = []
    for i in xrange(ntraj):
        therm_i = Thermostat({"nu_therm":params["nu_therm"], 
                              "NHC_size":params["NHC_size"], 
                              "Temperature":params["Temperature"], 
                              "thermostat_type":params["thermostat_type"]
                             })
        therm_i.set_Nf_t(nnucl)
        therm_i.set_Nf_r(0)
        therm_i.init_nhc()
        therm.append(therm_i)

    if verbose==1:
        print "therm =",therm 

    # Returned values:
    # therm - the list of thermostat objects
    return therm

