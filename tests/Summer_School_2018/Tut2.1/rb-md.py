#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

###################################################################
# Tutorial: Here is another version of MD: RB-MD via Systems objects
# Now, lets make it in a modular way
# and also add cooling
###################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


kb = (1.9872065E-3 / 627.5094709)  # in Ha now

##############################################################

def init_md(syst, mol, el, ham):
    syst.set_atomic_q(mol.q)  # mol -> syst

    E_kin = 0.0
    for n in xrange(syst.Number_of_fragments):
        E_kin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    syst.zero_forces_and_torques()

    # Compute forces and 
    syst.extract_atomic_q(mol.q) # syst -> mol
    E_pot = compute_potential_energy(mol, el, ham, 1)  # 1 - FSSH forces
    compute_forces(mol, el, ham, 1)
    syst.set_atomic_f(mol.f)  # mol -> syst
    syst.update_fragment_forces_and_torques();

    E_tot = E_pot + E_kin;

    return [E_kin, E_pot, E_tot]


def md_step(syst, mol, el, ham, dt, integrator):

    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    ps = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Propagate translational DOFs 
        syst.Fragments[n].Group_RB.shift_position(dt * syst.Fragments[n].Group_RB.rb_p * syst.Fragments[n].Group_RB.rb_iM);
 
        # Propagate rotational DOFs
        if integrator=="Jacobi":
            syst.Fragments[n].Group_RB.propagate_exact_rb(dt)
        elif integrator=="DLML":
            ps = syst.Fragments[n].Group_RB.propagate_dlml(dt)

    # Update atomic coordinates
    for n in xrange(syst.Number_of_fragments):
        syst.update_atoms_for_fragment(n)

    #=========== Potential ===================
    syst.zero_forces_and_torques()
    syst.extract_atomic_q(mol.q)

    epot = compute_potential_energy(mol, el, ham, 1)  # 1 - FSSH forces
    compute_forces(mol, el, ham, 1)

    syst.set_atomic_f(mol.f)
    syst.update_fragment_forces_and_torques()


    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    etot = ekin+epot
    return [ekin, epot, etot]




def run_rb_md(params, rnd):
    ##############################################################

    # Create Universe and populate it
    U = Universe()
    verbose = 0
    LoadPT.Load_PT(U, params["elements_file"], verbose)

    # Create force field
    uff = ForceField(params["ff"])
    LoadUFF.Load_UFF(uff)
    verb, assign_rings = 0, 1

    # Create system: load from a file and analyze
    syst = System()
    LoadMolecule.Load_Molecule(U, syst, os.getcwd()+params["input_structure"], "pdb")    
    syst.determine_functional_groups(1)  # 
    syst.show_atoms()
    syst.init_fragments()
    nat = syst.Number_of_atoms;   print "Number of atoms in the system = ", nat
    atlst1 = range(1,syst.Number_of_atoms+1)

    # Creating Hamiltonian
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)  
    ham.show_interactions_statistics()

    # Associate the system with the Hamiltonian and perform first calculaitons
    ham.set_system(syst)
    ham.compute()
    print "Initial point:"
    print "Energy = ", ham.H(0,0), " a.u."
    print "Gradient on coordinate 1 (x of first atom) = ", ham.dHdq(0,0,0), " a.u."
    print "Gradient on coordinate 3 (x of second atom) = ", ham.dHdq(0,0,3), " a.u."



    #== Prepare to run molecular dynamics ==
    # Electornic DOFs
    el = Electronic(1,0)

    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    print "q= ", mol.q[0]
    print "p= ", mol.p[0]
    print "f= ", mol.f[0]
    print "mass= ", mol.mass[0]
    
    init_md(syst, mol, el, ham)
    syst.init_fragments()
    syst.show_fragments()


    #========================= Cooling =================================
    f = open("_en_cooling.txt","w")
    f.close()

    for i in xrange(params["ncycles_cooling"]):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_cooling.xyz",i)

        for j in xrange(params["nsteps_cooling"]):
            ekin, epot, etot = md_step(syst, mol, el, ham, params["dt_cooling"], params["integrator"])
            Tcurr = 2.0 * ekin / (3.0 * nat * kb)

        syst.cool()

        f = open("_en_cooling.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  Temperature= %8.5f \n" % (i, ekin, epot, etot, Tcurr))
        f.close()

    #========================= Production MD =================================
    f = open("_en_md.txt","w")
    f.close()

    syst.init_fragment_velocities(params["T"],params["Tot_P"], params["Tot_L"], rnd)

    for i in xrange(params["ncycles_md"]):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)

        for j in xrange(params["nsteps_md"]):
            ekin, epot, etot = md_step(syst, mol, el, ham, params["dt_md"], params["integrator"])
            Tcurr = 2.0 * ekin / (3.0 * nat * kb)

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  Temperature= %8.5f \n" % (i, ekin, epot, etot, Tcurr))
        f.close()





########################## Simulations #################################

rnd = Random()

params = {}
params["elements_file"] = "elements.dat"
params["ff"] = { "bond_functional":"Harmonic",  "angle_functional":"Fourier", "dihedral_functional":"General0",
                 "oop_functional":"Fourier",   "mb_functional":"LJ_Coulomb","R_vdw_on":10.0,"R_vdw_off":15.0 }

params["input_structure"] = "/23waters.ent"
#params["input_structure"] = "/23waters-aa.ent"

params["integrator"] = "DLML"

params["dt_cooling"] = 20.0
params["ncycles_cooling"] = 100
params["nsteps_cooling"] = 50

params["T"] = 200.0
params["Tot_P"] = VECTOR(0.0, 0.0, 0.0)
params["Tot_L"] = VECTOR(0.0, 0.0, 0.0)
params["dt_md"] = 20.0
params["ncycles_md"] = 100
params["nsteps_md"] = 100

run_rb_md(params, rnd)
