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




##############################################################

# Create Universe and populate it
U = Universe()
verbose = 0
LoadPT.Load_PT(U, "elements.dat", verbose)


# Create force field
uff = ForceField({"bond_functional":"Harmonic",
                  "angle_functional":"Fourier",
                  "dihedral_functional":"General0",
                  "oop_functional":"Fourier",
                  "mb_functional":"LJ_Coulomb","R_vdw_on":10.0,"R_vdw_off":15.0 })



LoadUFF.Load_UFF(uff)
verb = 0
assign_rings = 1



#======= System ==============
#for i in range(1,13):
for i in [1]:
    print "=================== System ",i,"======================="

    syst = System()
    LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/Clusters/23waters.ent", "pdb")
    #Load_Molecule(U, syst, os.getcwd()+"/Clusters/23waters_noq.ent", "pdb")

    
    syst.determine_functional_groups(1)  # 
    syst.show_atoms()

    syst.init_fragments()
#    syst.show_fragments()
#    syst.show_molecules()

#    syst.
    

    #syst.print_xyz("molecule1.xyz",1)

    print "Number of atoms in the system = ", syst.Number_of_atoms
    atlst1 = range(1,syst.Number_of_atoms+1)


    # Creating Hamiltonian
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)  

    ham.show_interactions_statistics()

   
    ham.set_system(syst)
    ham.compute()



#    sys.exit(0)

    print "Energy = ", ham.H(0,0), " a.u."
    print "Force 1 = ", ham.dHdq(0,0,0), " a.u."
    print "Force 3 = ", ham.dHdq(0,0,3), " a.u."


#    sys.exit(0)

    #--------------------- Molecular dynamics ----------------------
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


#    sys.exit(0)



    #=================== Propagation ====================

    integrator = "DLML"

    f = open("_en_traj.txt","w")
    dt = 20.0 # = 0.5 fs
    for i in xrange(100):

        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_traj.xyz",i)

        for j in xrange(100):

            ekin = 0.0
            for n in xrange(syst.Number_of_fragments):
                # Linear momentum propagation:
                #syst.Fragments[n].Group_RB.scale_linear_(0.5*dt)
                syst.Fragments[n].Group_RB.apply_force(0.5*dt)

                # Angular momentum propagation:
                #syst.Fragments[n].Group_RB.scale_angular_(0.5*dt)
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
                #syst.Fragments[n].Group_RB.scale_linear_(0.5*dt)
                syst.Fragments[n].Group_RB.apply_force(0.5*dt)

                # Angular momentum propagation:
                #syst.Fragments[n].Group_RB.scale_angular_(0.5*dt)
                syst.Fragments[n].Group_RB.apply_torque(0.5*dt)
                
                ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )


        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, ekin, epot, ekin+epot))

        


      
    f.close()



