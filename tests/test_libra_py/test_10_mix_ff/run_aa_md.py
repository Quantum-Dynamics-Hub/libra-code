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

###################################################################
# Starting with the file from the test_hamiltonian_mm/test_mm6a.py
#
###################################################################

import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

rnd = Random()

def main():

    #--------------------- Initialization ----------------------

    # Create Universe and populate it
    U = Universe(); LoadPT.Load_PT(U, "elements.dat")

    # Create force field
    uff = ForceField({"bond_functional":"Harmonic"}) #, "mb_functional":"LJ_Coulomb","R_vdw_on":10.0,"R_vdw_off":15.0 })
    LoadUFF.Load_UFF(uff,"uff.dat")


    Xff = ForceField({"angle_functional":"Fourier"})  # no FF name, so we'll rely on the manually-defined atom type    


    class atomrecord:
        pass

    ff = atomrecord()
    ff.Atom_ff_int_type = 1
    ff.Atom_ff_type = "H_"
    atom_record = Atom_Record()
    atom_record.set(ff)
    Xff.Add_Atom_Record(atom_record)


    ff = atomrecord()
    ff.Atom_ff_int_type = 2
    ff.Atom_ff_type = "H_b"
    atom_record = Atom_Record()
    atom_record.set(ff)
    Xff.Add_Atom_Record(atom_record)


    class anglerecord:
        pass

    ff = anglerecord()
    ff.Atom1_ff_type    = "H_"
    ff.Atom2_ff_type    = "H_b"
    ff.Atom3_ff_type    = "H_b"
    ff.Angle_k_angle    = 100.0
    ff.Angle_theta_eq   = 120.0
    angle_record = Angle_Record()
    angle_record.set(ff)
    angle_record.show_info()
    Xff.Add_Angle_Record(angle_record)

    ff = anglerecord()
    ff.Atom1_ff_type    = "H_b"
    ff.Atom2_ff_type    = "H_b"
    ff.Atom3_ff_type    = "H_b"
    ff.Angle_k_angle    = 100.0
    ff.Angle_theta_eq   = 120.0
    angle_record = Angle_Record()
    angle_record.set(ff)
    angle_record.show_info()
    Xff.Add_Angle_Record(angle_record)

    Xff.show_atom_records()
    Xff.show_angle_records()




    # Create molecular system and initialize the properties
    syst = System()
    LoadMolecule.Load_Molecule(U, syst, "Hn.ent", "pdb")
    syst.determine_functional_groups(0)  # do not assign rings
    syst.init_fragments()
    print "Number of atoms in the system = ", syst.Number_of_atoms
    atlst1 = range(1,syst.Number_of_atoms+1)


    # Creating Hamiltonian and initialize it
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, 1, 0)  # 0 - verb, 0 - assign_rings
    ham.show_interactions_statistics()
    ham.show_interactions()


    #======= Now here, we'll add extra interactions due another FF ============
    # keep in mind that the FF atom types in syst are already set to those of UFF
    # so, we don't need to define extra types
    #
    # Also, we need to modify the coordination of all atoms
    for i in xrange(syst.Number_of_atoms):
        syst.Atoms[i].Atom_coordination = 1
        syst.Atoms[i].is_Atom_coordination = 1

    ham.set_interactions_for_atoms(syst, atlst1, atlst1, Xff, 1, 0)  # 0 - verb, 0 - assign_rings    
    ham.show_interactions_statistics() 
    ham.show_interactions()


    # Bind Hamiltonian and the system   
    ham.set_system(syst);   ham.compute();   print "Energy = ", ham.H(0,0), " a.u."

    # Electronic DOFs
    el = Electronic(1,0)

    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)

    # Initialize MD variables
    nve_md.nve_md_init(syst, mol, el, ham)

    #=================== Propagation ====================
    integrator = "DLML"

    ########################## Production MD #################################

    syst.init_atom_velocities(300.0, rnd)

    f = open("_en_md.txt","w")
    f.close()
    dt = 20.0 

    for i in xrange(100):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)

        for j in xrange(25):
            ekin, epot, etot = nve_md.nve_md_step(syst, mol, el, ham, dt, integrator)

        syst.cool()

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, ekin, epot, etot))
        f.close()


main()

