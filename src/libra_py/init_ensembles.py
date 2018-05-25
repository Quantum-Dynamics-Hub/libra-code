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

import LoadMolecule
import LoadPT




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


def init_systems(ntraj, U, filename, format, verbose=0):
    ##
    # This function creates a list of System objects and loads the content into them
    # 
    # \param[in] ntraj The number of trajectories in the ensemble
    # \param[in] U Is the Universe object 
    # \param[in] filename The name of the file containing molecular structure
    # \param[in] A string defining the file format (see LoadMolecule for more details)
    # \param[in] verbose The parameter controlling the amount of debug printing: 0 (default)
    # - no printing, 1 - print info

    syst = []
    for tr in xrange(ntraj):
        syst.append( System() )
        LoadMolecule.Load_Molecule(U, syst[tr], filename, format)
        syst[tr].determine_functional_groups(0)  # do not assign rings
        syst[tr].init_fragments()

        if tr==0 and verbose==1:
            print "Number of atoms in the system = ", syst[tr].Number_of_atoms
            print "Number of bonds in the system = ", syst[tr].Number_of_bonds
            print "Number of angles in the system = ", syst[tr].Number_of_angles
            print "Number of dihedrals in the system = ", syst[tr].Number_of_dihedrals
            print "Number of impropers in the system = ", syst[tr].Number_of_impropers

    return syst


def init_systems2(ntraj, filename, format, rnd, T, sigma, data_file="elements.dat", verbose=0):
    ##
    # This function creates a list of System objects and loads the content into them, it also randomizes
    # the original input from the coordinate file (by adding coordinate displacements and setting velocities)
    # 
    # \param[in] ntraj The number of trajectories in the ensemble
    # \param[in] filename The name of the file containing molecular structure
    # \param[in] format A string defining the file format (see LoadMolecule for more details)
    # \param[in] rnd   Random number generator object
    # \param[in] T     Target temperature used to initialize momenta of atoms. units: K
    # \param[in] sigma The magnitude of a random displacement of each atom from its center, units: Bohr
    # \param[in] data_file The file containing information about elements.
    # \param[in] verbose The parameter controlling the amount of debug printing: 0 (default)
    # - no printing, 1 - print info


    # Create Universe and populate it
    U = Universe();   LoadPT.Load_PT(U, data_file, 0)
    syst = []

    for tr in xrange(ntraj):
        syst.append( System() )
        LoadMolecule.Load_Molecule(U, syst[tr], filename, format)

        for at in xrange(syst[tr].Number_of_atoms):
            # Random shift
            shift = VECTOR(sigma*rnd.normal(), sigma*rnd.normal(), sigma*rnd.normal() )

            # Apply it
            syst[tr].Atoms[at].Atom_RB.shift_position(shift)


        # Initialize random velocity at T(K) using normal distribution
        syst[tr].init_atom_velocities(T,rnd)


        syst[tr].determine_functional_groups(0)  # do not assign rings
        syst[tr].init_fragments()

        if tr==0 and verbose==1:
            print "Number of atoms in the system = ", syst[tr].Number_of_atoms
            print "Number of bonds in the system = ", syst[tr].Number_of_bonds
            print "Number of angles in the system = ", syst[tr].Number_of_angles
            print "Number of dihedrals in the system = ", syst[tr].Number_of_dihedrals
            print "Number of impropers in the system = ", syst[tr].Number_of_impropers

    return syst


                      
def init_mm_Hamiltonians(syst, ff, verbose=0):
    ##
    # This function creates a list of MM Hamiltonians and bind them to the corresponding systems
    # 
    # \param[in] syst The list of System objects 
    # \param[in] ff The Force Field object
    # \param[in] verbose The parameter controlling the amount of debug printing: 0 (default)
    # - no printing, 1 - print info

    ntraj = len(syst)
    ham = []

    # Creating Hamiltonian and initialize it
    for tr in xrange(ntraj):
        ham.append( Hamiltonian_Atomistic(1, 3*syst[tr].Number_of_atoms) )
        ham[tr].set_Hamiltonian_type("MM")


        atlst1 = range(1,syst[tr].Number_of_atoms+1)
        ham[tr].set_interactions_for_atoms(syst[tr], atlst1, atlst1, ff, 1, 0)  # 0 - verb, 0 - assign_rings
           
        # Bind MM-Hamiltonian and the system   
        ham[tr].set_system(syst[tr]);
        ham[tr].compute(); 

        if tr==0 and verbose==1:
            print "Printing the Hamiltonian info for the first trajectory"
            ham[tr].show_interactions_statistics()
            print "Energy = ", ham.H(0,0), " a.u."

    return ham



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

