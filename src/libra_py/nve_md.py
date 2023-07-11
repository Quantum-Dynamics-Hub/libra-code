#*********************************************************************************
#* Copyright (C) 2016-2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: nve_md
   :platform: Unix, Windows
   :synopsis: This module implements a number of auxiliary functions for simpler
       implementation of classical MD simulations
.. moduleauthor:: Alexey V. Akimov

"""

import sys
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
from . import units
from . import LoadPT
from . import LoadUFF


##############################################################

def syst2xyz(syst):
    """

    This function generates the xyz string representing the system

    Args:
        syst ( System object ): the chemical system
 
    Returns:
        string:  res:  the coordinates of the system in an xyz format

    """
     
    nat = syst.Number_of_atoms
    res = " %i \n\n" % (nat)
    for i in range(0,nat):
        x = syst.Atoms[i].Atom_RB.rb_cm.x
        y = syst.Atoms[i].Atom_RB.rb_cm.y
        z = syst.Atoms[i].Atom_RB.rb_cm.z
        res = res + "%s  %5.3f %5.3f %5.3f\n" % (syst.Atoms[i].Atom_element, x, y, z)

    return res    


def nve_md_init(syst, mol, el, ham):
    """

    This function initializes all the necessary variables for classical MD

    Args:
    syst ( System object ): The object containing all the information about the chemical system.
        This object is updated
    mol ( Nuclear object ): The nuclear DOF. This object is updated
    el ( Electronic object ): The Electronic object describing which PES to use (in classical MD
        it is typically the ground state)
    ham ( Hamiltonian object ): The Hamiltonian object representing all the interactions

    Returns: 
        (double, double, double): E_kin, E_pot, E_tot:

            * E_kin: kinetic energy [ units: Ha ]
            * E_pot: potential energy [ units: Ha ]
            * E_tot: total energy [ units: Ha ]

        This function also updates the variables of the input variable ```syst```
     
    Note:

        Operations:
        1. syst -> mol 
        2. compute kinetic energy
        3. compute forces and potential energy
        4. compute total energy
        5. initialize fragment variables, if any

    """

    # 1. syst-> mol
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    # 2. compute kinetic energy
    E_kin = 0.0
    for n in range(0,syst.Number_of_fragments):
        E_kin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    # 3. Compute forces and potential energy
    syst.zero_forces_and_torques()
    E_pot = compute_forces(mol, el, ham, 1)  # - FSSH forces
    syst.set_atomic_f(mol.f)  # mol -> syst
    syst.update_fragment_forces_and_torques();

    # 4. Compute total energy
    E_tot = E_pot + E_kin 

    # 5. Initialize fragment variables
    syst.init_fragments()

    return E_kin, E_pot, E_tot




def nve_md_step(syst, mol, el, ham, params): 
    """
    This function performs a classical NVE MD step

    Args: 
        syst ( System object ): The object containing all the information about the chemical system.
            This object is updated
        mol ( Nuclear object ): The nuclear DOF. This object is updated
        el ( Electronic object ): The Electronic object describing which PES to use (in classical MD
            it is typically the ground state)
        ham ( Hamiltonian object ): The Hamiltonian object representing all the interactions
        params ( dictionary ): The parameters controlling the execution of the dynamics
        
            * **params["dt"]** ( double ): integration timestep [ units: a.u., default: 20.0 ]
            * **params["integrator"]** ( string ): The rigid-body MD integrator [ default: "DLML" ]
            * **params["fixed_fragment_translation"]** ( list of ints ): the indices (starting from 0)
                of the fragments whose translational DOFs are frozen [ default: empty ]
            * **params["fixed_fragment_rotation"]** ( list of ints ): the indices (starting from 0)
                of the fragments whose rotational DOFs are frozen [ default: empty ]


    Returns: 
        (double, double, double): E_kin, E_pot, E_tot:

            * E_kin: kinetic energy [ units: Ha ]
            * E_pot: potential energy [ units: Ha ]
            * E_tot: total energy [ units: Ha ]

        This function also updates the variables of the input variable ```syst```
     
    Note:
        Operations:
        1. propagate rotational and translational momenta for 0.5 of dt
        2. propagate rotational and translational coordinates of all fragments for dt
        3. update atomic positions
        4. compute potential energy and forces on atoms
        5. update fragmental forces and torques
        6. propagate rotational and translational momenta for another 0.5 of dt

    """

    #=== Setup parameters =======
    critical_params = [ ] 
    default_params = {"dt":20.0, "integrator":"DLML", "fixed_fragment_translation":[], 
                      "fixed_fragment_rotation":[]
                     }
    comn.check_input(params, default_params, critical_params)

    dt = params["dt"]
    integrator = params["integrator"]
    fixed_tr = params["fixed_fragment_translation"]
    fixed_rot = params["fixed_fragment_rotation"]


    # 1. propagate rotational and translational DOF for 0.5 of dt
    ekin = 0.0
    for n in range(0,syst.Number_of_fragments):
        # Linear momentum propagation:
        if n not in fixed_tr:
            syst.Fragments[n].Group_RB.apply_force(0.5*dt)

        # Angular momentum propagation:
        if n not in fixed_rot:
            syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        

        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    # 2. propagate rotational and translational coordinates of all fragments for dt
    ps = 0.0
    for n in range(0,syst.Number_of_fragments):
        # Propagate translational DOFs 
        if n not in fixed_tr:
            syst.Fragments[n].Group_RB.shift_position(dt * syst.Fragments[n].Group_RB.rb_p * syst.Fragments[n].Group_RB.rb_iM);
 
        # Propagate rotational DOFs
        if n not in fixed_rot:
            if integrator=="Jacobi":
                syst.Fragments[n].Group_RB.propagate_exact_rb(dt)
            elif integrator=="DLML":
                ps = syst.Fragments[n].Group_RB.propagate_dlml(dt)

    # 3. update atomic positions
    for n in range(0,syst.Number_of_fragments):
        syst.update_atoms_for_fragment(n)

    # 4. compute potential energy and forces on atoms
    syst.zero_forces_and_torques()
    syst.extract_atomic_q(mol.q)  # syst -> mol
    epot = compute_forces(mol, el, ham, 1)  # 1 - FSSH forces
    syst.set_atomic_f(mol.f)   # mol -> syst

    # 5. update fragmental forces and torques
    syst.update_fragment_forces_and_torques()

    # 6. propagate rotational and translational momenta for another 0.5 of dt
    ekin = 0.0
    for n in range(0,syst.Number_of_fragments):
        # Linear momentum propagation:
        if n not in fixed_tr:
            syst.Fragments[n].Group_RB.apply_force(0.5*dt)

        # Angular momentum propagation:
        if n not in fixed_rot:
            syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    etot = ekin+epot
    return ekin, epot, etot


#syst = System()
#LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/Clusters/23waters.ent", "pdb")    

def optimize_syst(syst, params):
    """
    A function to optimize the geometry of the system 

    Args:
        syst ( System object ): represents the chemical object
        params ( dictionary ): control parameters:

            * **params["anneal_schedule"]** ( list ): the annealing schedule. Each element of the list
                consist of 3 elements: dt, ncycles, nsteps, where:
                dt - the timestep for integration [ units: a.u. ]
                ncycles - the number of cycles of annealing with dt 
                Each annealing cycle consists of ```nsteps```  steps of NVE MD steps followed by the cooling
                Cooling just resets all the momenta and angular momenta to zero
                Example:
                    params["anneal_schedule"] = [ [1.0, 100, 10], [20.0, 100, 100] ] means:
                    First do 100 cycles of annealing: 10 MD steps with dt = 1 a.u. each followed by cooling
                    Second do 100 cycles of annealing: 100 MD steps with dt = 20 a.u. each followed by cooling                
 
            * **params["elements_file"]** ( string ): The file that contains properties of the elements 
                it is needed for construction of the Universe [ default: "elements.dat" ]
            * **params["cooling_out1"]** ( Boolean ): Whether to print out the energies along the simulated
                cooling protocol. If selected, the info is printed out to the file "_en_cooling.txt" [ default: False ]
            * **params["cooling_out2"]** ( Boolean ): Whether to print out the coordinates in the xyz format.
                If selected, the info is printed out to the "_mol_cooling.xyz" file [ default: False]
            [ default: False ]


            SeeAlso: is ```nve_md_step``` 
            * **params["dt"]** ( double ): integration timestep [ units: a.u., default: 20.0 ]
            * **params["integrator"]** ( string ): The rigid-body MD integrator [ default: "DLML" ]
            * **params["fixed_fragment_translation"]** ( list of ints ): the indices (starting from 0)
                of the fragments whose translational DOFs are frozen [ default: empty ]
            * **params["fixed_fragment_rotation"]** ( list of ints ): the indices (starting from 0)
                of the fragments whose rotational DOFs are frozen [ default: empty ]



    """

    critical_params = [ ] 
    default_params = {"elements_file":"elements.dat", "cooling_out1":False, "cooling_out2":False,
                      "anneal_schedule":[[20.0, 100, 10]],
                      "fixed_fragment_translation":[], 
                      "fixed_fragment_rotation":[]
                     }
    comn.check_input(params, default_params, critical_params)
    
    elements_file = params["elements_file"]
    cooling_out1 = params["cooling_out1"]
    cooling_out2 = params["cooling_out2"]
    anneal_schedule = params["anneal_schedule"]


    # Create Universe and populate it
    U = Universe()
    verbose = 0
    LoadPT.Load_PT(U, elements_file, verbose)

    # Create force field
    uff = ForceField({"bond_functional":"Harmonic",
                      "angle_functional":"Fourier",
                      "dihedral_functional":"General0",
                      "oop_functional":"Fourier",
                      "mb_functional":"LJ_Coulomb","R_vdw_on":10.0,"R_vdw_off":15.0
                     })
    LoadUFF.Load_UFF(uff)

    syst.determine_functional_groups(1)  # 
    syst.init_fragments()
    atlst1 = list(range(1,syst.Number_of_atoms+1))


    # Creating Hamiltonian
    verb, assign_rings = 0, 1
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)     
    ham.set_system(syst)
    ham.compute()

    # Electronic DOF
    el = Electronic(1,0)

    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)  # syst -> mol
    syst.extract_atomic_p(mol.p)  # syst -> mol
    syst.extract_atomic_f(mol.f)  # syst -> mol
    syst.extract_atomic_mass(mol.mass)

    #init_md(syst, mol, el, ham)
    nve_md_init(syst, mol, el, ham)
    syst.init_fragments()


    ########################## Cooling #################################
    params1 = dict(params)

    if cooling_out1:
        f = open("_en_cooling.txt","w")
        f.close()

    for anneal_item in anneal_schedule:
        params1["dt"] = anneal_item[0]
        ncycles = anneal_item[1]
        nsteps = anneal_item[2]

        for i in range(0,ncycles):  # the number of cycles 
            if cooling_out2:
                syst.set_atomic_q(mol.q)  # mol -> syst - probably not needed
                syst.print_xyz("_mol_cooling.xyz",i)

            for j in range(0,nsteps):                
                ekin, epot, etot = nve_md_step(syst, mol, el, ham, params1)  # the number of steps 


            syst.cool()

            if cooling_out1:
                f = open("_en_cooling.txt","a")
                f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, ekin, epot, etot))
                f.close()


