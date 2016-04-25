#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file md.py 
# This module implements a number of auxiliary functions for simpler implementation of 
# classical MD simulations

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


##############################################################

def nve_md_init(syst, mol, el, ham):
## This function initializes all the necessary variables
# \param[in,out] syst The chemsystem object containing all the information about the system
# \param[out] mol The Nuclear object containing nuclear DOF
# \param[in] el The Electronic object describing which PES to use (in classical MD it is typically the ground state)
# \param[in] ham The Hamiltonian object representing all the interactions
#
# Operations:
# 1. syst -> mol 
# 2. compute kinetic energy
# 3. compute forces and potential energy
# 4. compute total energy
# 5. initialize fragment variables
#
# Returns: E_kin, E_pot, E_tot
# where:
# E_kin - kinetic energy
# E_pot - potential energy
# E_tot - total energy

    # 1. syst-> mol
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    # 2. compute kinetic energy
    E_kin = 0.0
    for n in xrange(syst.Number_of_fragments):
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




def nve_md_step(syst, mol, el, ham, dt=20.0, integrator="DLML"):
## This function performs a classical NVE MD step
# \param[in,out] syst The chemsystem object containing all the information about the system
# \param[in,out] mol The Nuclear object containing nuclear DOF
# \param[in] el The Electronic object describing which PES to use (in classical MD it is typically the ground state)
# \param[in] ham The Hamiltonian object representing all the interactions
# \param[in] dt Time step, in atomic units
# \param[in] integrator The rigid-body MD integrator
#
# Operations:
# 1. propagate rotational and translational momenta for 0.5 of dt
# 2. propagate rotational and translational coordinates of all fragments for dt
# 3. update atomic positions
# 4. compute potential energy and forces on atoms
# 5. update fragmental forces and torques
# 6. propagate rotational and translational momenta for another 0.5 of dt
#
# Returns: E_kin, E_pot, E_tot
# where:
# E_kin - kinetic energy
# E_pot - potential energy
# E_tot - total energy


    # 1. propagate rotational and translational DOF for 0.5 of dt
    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    # 2. propagate rotational and translational coordinates of all fragments for dt
    ps = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Propagate translational DOFs 
        syst.Fragments[n].Group_RB.shift_position(dt * syst.Fragments[n].Group_RB.rb_p * syst.Fragments[n].Group_RB.rb_iM);
 
        # Propagate rotational DOFs
        if integrator=="Jacobi":
            syst.Fragments[n].Group_RB.propagate_exact_rb(dt)
        elif integrator=="DLML":
            ps = syst.Fragments[n].Group_RB.propagate_dlml(dt)

    # 3. update atomic positions
    for n in xrange(syst.Number_of_fragments):
        syst.update_atoms_for_fragment(n)

    # 4. compute potential energy and forces on atoms
    syst.zero_forces_and_torques()
    syst.extract_atomic_q(mol.q)
    epot = compute_forces(mol, el, ham, 1)  # 1 - FSSH forces
    syst.set_atomic_f(mol.f)

    # 5. update fragmental forces and torques
    syst.update_fragment_forces_and_torques()

    # 6. propagate rotational and translational momenta for another 0.5 of dt
    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    etot = ekin+epot
    return ekin, epot, etot


