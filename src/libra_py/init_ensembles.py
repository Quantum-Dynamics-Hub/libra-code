#*********************************************************************************
#* Copyright (C) 2016-2019 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: init_ensembles
   :platform: Unix, Windows
   :synopsis: This module implements functions that allocate memory for ensembles of trajectories
.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import LoadMolecule
from . import LoadPT


def init_ext_hamiltonians(ntraj, nnucl, nel, verbose=0):
    """

    This function allocates memory for matrices (external objects). The function
    also creates an array of External_Hamiltonian objects and bind them to the external matrices

    Args:
        ntraj ( int ): The number of trajectories in the ensemble
        nnucl ( int ): The number of nuclear DOF ( usually, = 3* Number_of_atoms)
        nel ( int ): The number of electronic DOF ( usually the number of electronic states)
        verbose ( int ): The parameter controlling the amount of debug printing: 

            - 0: no printing [ default ]
            - 1: do print info
    Returns:
        tuple: (ham, ham_adi, d1ham_adi, ham_vib), where:

            * ham ( list of ntraj External_Hamiltonian objects ): "External Hamiltonian" objects,
                one per trajectory. These are the "controller" objects that organize the flow of
                computations and store the results in the matrices bound to them.
            * ham_adi ( list of ntraj MATRIX(nel, nel) ): Electronic Hamiltonian matrices 
                coupled to the ham objects
            * d1ham_adi ( list ntraj lists of nnucl MATRIX(nel, nel) objects ): derivatives of 
                electronic Hamiltonians w.r.t. nuclear DOFs for each trajectory. Also coupled to 
                the ham objects for the corresponding trajectories.
            * ham_vib ( list of ntraj CMATRIX(nel, nel) ): Vibronic Hamiltonian matrices 
                coupled to the ham objects
     
    """

    # Create an list of External Hamiltonian objects
    ham = []
    for i in range(0,ntraj):
        ham_i = Hamiltonian_Extern(nel, nnucl) 
        ham_i.set_rep(1)            # adiabatic representation
        ham_i.set_adiabatic_opt(0)  # use the externally-computed 
                                    # adiabatic electronic Hamiltonian and derivatives
        ham_i.set_vibronic_opt(0)   # use the externally-computed vibronic Hamiltonian and derivatives
        ham.append(ham_i)

    # bind actual matrices to the External_Hamiltonian objects
    ham_adi = []; d1ham_adi = []; ham_vib = []


    for i in range(0,ntraj):
        # create external matrix for adiabatic Hamiltonian
        ham_adi_i = MATRIX(nel, nel)
        ham_adi.append(ham_adi_i)
        # bind adiabatic Hamiltonian
        ham[i].bind_ham_adi(ham_adi[i])

        # create external matrices for the derivatives
        d1ham_adi_i= MATRIXList()
        for k in range(0,nnucl):
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
        print("Addresses of the working matrices")
        for i in range(0,ntraj):
            print("%i th Hamiltonian is" % i)
            print("ham_adi = ", ham_adi[i])
            print("d1ham_adi = ", d1ham_adi[i])
            for k in range(0,nnucl):
                print("d1ham_adi[",k,"]= ", d1ham_adi[i][k])
            print("ham_vib = ", ham_vib[i])

    return ham, ham_adi, d1ham_adi, ham_vib


def init_systems(ntraj, U, filename, format, verbose=0):
    """
    This function creates a list of System objects and loads the content into them

    Args:
        ntraj ( int ): The number of trajectories in the ensemble
        U ( Universe object ): The container of the fundamental parameters of atoms
        filename ( string ): The file containg the atomistic system 
        format ( string ): The format of the file to load (see LoadMolecule for possible formats)
        verbose ( int ): The parameter controlling the amount of debug printing: 

            - 0: no printing [ default ]
            - 1: do print info

    Returns:
        (list of ntraj System objects): the ensemble of the objects defining the replicas of the system

    """

    syst = []
    for tr in range(0,ntraj):
        syst.append( System() )
        LoadMolecule.Load_Molecule(U, syst[tr], filename, format)
        syst[tr].determine_functional_groups(0)  # do not assign rings
        syst[tr].init_fragments()

        if tr==0 and verbose==1:
            print("Number of atoms in the system = ", syst[tr].Number_of_atoms)
            print("Number of bonds in the system = ", syst[tr].Number_of_bonds)
            print("Number of angles in the system = ", syst[tr].Number_of_angles)
            print("Number of dihedrals in the system = ", syst[tr].Number_of_dihedrals)
            print("Number of impropers in the system = ", syst[tr].Number_of_impropers)

    return syst


def init_systems2(ntraj, filename, format, rnd, T, sigma, data_file="elements.dat", verbose=0):
    """
    This function creates a list of System objects and loads the content into them, it also randomizes
    the original input from the coordinate file (by adding coordinate displacements and setting velocities)

    Args:
        ntraj ( int ): The number of trajectories in the ensemble
        U ( Universe object ): The container of the fundamental parameters of atoms
        filename ( string ): The file containg the atomistic system 
        format ( string ): The format of the file to load (see LoadMolecule for possible formats)
        rnd ( Random object ): random number generator object
        T ( double ): Target temperature used to initialize momenta of atoms. [ units: K ]
        sigma ( double ): The magnitude of a random displacement of each atom from its center [ units: Bohr ]
        data_file ( double ): The file containing information about elements [ default: "elements.dat" ]
        verbose ( int ): The parameter controlling the amount of debug printing: 

            - 0: no printing [ default ]
            - 1: do print info

    Returns:
        (list of ntraj System objects): the ensemble of the objects defining the replicas of the system

    """

    # Create Universe and populate it
    U = Universe();   LoadPT.Load_PT(U, data_file, 0)
    syst = []

    for tr in range(0, ntraj):
        syst.append( System() )
        LoadMolecule.Load_Molecule(U, syst[tr], filename, format)

        for at in range(0, syst[tr].Number_of_atoms):
            # Random shift
            shift = VECTOR(sigma*rnd.normal(), sigma*rnd.normal(), sigma*rnd.normal() )

            # Apply it
            syst[tr].Atoms[at].Atom_RB.shift_position(shift)


        # Initialize random velocity at T(K) using normal distribution
        syst[tr].init_atom_velocities(T,rnd)


        syst[tr].determine_functional_groups(0)  # do not assign rings
        syst[tr].init_fragments()

        if tr==0 and verbose==1:
            print("Number of atoms in the system = ", syst[tr].Number_of_atoms)
            print("Number of bonds in the system = ", syst[tr].Number_of_bonds)
            print("Number of angles in the system = ", syst[tr].Number_of_angles)
            print("Number of dihedrals in the system = ", syst[tr].Number_of_dihedrals)
            print("Number of impropers in the system = ", syst[tr].Number_of_impropers)

    return syst


                      
def init_mm_Hamiltonians(syst, ff, verbose=0):
    """

    This function creates a list of MM Hamiltonians and bind them to the corresponding systems
    
    Args:
        syst ( list of ntraj System objects ): self-explanatory
        ff ( ForceField object ): defines which force fields to construct for all systems
        verbose ( int ): The parameter controlling the amount of debug printing: 

            - 0: no printing [ default ]
            - 1: do print info

    Returns:
        (list of ntraj Hamiltonian_Atomistic objects): the ensemble of the objects defining
        the MM Hamiltonians that are also bound to the provided systems.

    """

    ntraj = len(syst)
    ham = []

    # Creating Hamiltonian and initialize it
    for tr in range(0,ntraj):
        ham.append( Hamiltonian_Atomistic(1, 3*syst[tr].Number_of_atoms) )
        ham[tr].set_Hamiltonian_type("MM")


        atlst1 = range(1,syst[tr].Number_of_atoms+1)
        ham[tr].set_interactions_for_atoms(syst[tr], atlst1, atlst1, ff, 1, 0)  # 0 - verb, 0 - assign_rings
           
        # Bind MM-Hamiltonian and the system   
        ham[tr].set_system(syst[tr]);
        ham[tr].compute(); 

        if tr==0 and verbose==1:
            print("Printing the Hamiltonian info for the first trajectory")
            ham[tr].show_interactions_statistics()
            print("Energy = ", ham.H(0,0), " a.u.")

    return ham




def init_mols(syst, ntraj, nnucl, verbose=0):
    """

    This function creates a list of Nuclear objects and initializes 
    them by extracting info from the provided systems - in "syst"

    Args:
        ntraj ( int ): The number of trajectories in the ensemble
        nnucl ( int ): The number of nuclear DOF ( usually, = 3* Number_of_atoms)
        verbose ( int ): The parameter controlling the amount of debug printing: 

            - 0: no printing [ default ]
            - 1: do print info
    Returns:
        (list of ntraj Nuclear objects): mol: self-explanatory

    """

    # Initialize nuclear variables
    mol = []
    for i in range(0,ntraj):
        mol.append(Nuclear(nnucl))
        syst[i].extract_atomic_q(mol[i].q)
        syst[i].extract_atomic_p(mol[i].p)
        syst[i].extract_atomic_f(mol[i].f)
        syst[i].extract_atomic_mass(mol[i].mass)

        if verbose==1:
            print(i,"mol=",mol[i])
            for k in range(0,syst[i].Number_of_atoms):
                print("mol[%d][%d]=%f,%f,%f"%(i,k,mol[i].q[3*k+0],mol[i].q[3*k+1],mol[i].q[3*k+2]) )

    return mol


def init_therms(ntraj, nnucl, params, verbose=0):
    """

    This function creates a list of Thermostat objects and initializes 
    them according to given parameters

    Args:
        ntraj ( int ): The number of trajectories in the ensemble
        nnucl ( int ): The number of nuclear DOF ( usually, = 3* Number_of_atoms)
        verbose ( int ): The parameter controlling the amount of debug printing: 

            - 0: no printing [ default ]
            - 1: do print info
    Returns:
        (list of ntraj Thermostat objects): therm: self-explanatory

    """

    # Initialize Thermostat object    
    therm = []
    for i in range(0,ntraj):
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
        print("therm =",therm)

    return therm

