#*********************************************************************************
#* Copyright (C) 2016-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

    

def main():

    rnd = Random()

    #--------------------- Initialization ----------------------

    # Create Universe and populate it
    U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.txt")

    # Create force field
    uff = ForceField({"bond_functional":"Harmonic",        "angle_functional":"Fourier",
                      "dihedral_functional":"General0",    "oop_functional":"Fourier",
                      "mb_functional":"LJ_Coulomb","R_vdw_on":40.0,"R_vdw_off":55.0 })
    LoadUFF.Load_UFF(uff,"uff.dat")


    # Create molecular system and initialize the properties
    syst = System()

    LoadMolecule.Load_Molecule(U, syst, "Pc-C60.ent", "pdb")

    syst.determine_functional_groups(0)  # do not assign rings
    syst.init_fragments()
    print "Number of atoms in the system = ", syst.Number_of_atoms
    print "Number of bonds in the system = ", syst.Number_of_bonds
    print "Number of angles in the system = ", syst.Number_of_angles
    print "Number of dihedrals in the system = ", syst.Number_of_dihedrals
    print "Number of impropers in the system = ", syst.Number_of_impropers
    atlst1 = range(1,syst.Number_of_atoms+1)

    # Creating Hamiltonian and initialize it
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, 1, 0)  # 0 - verb, 0 - assign_rings
    ham.show_interactions_statistics()


    # Bind Hamiltonian and the system   
    ham.set_system(syst);   ham.compute();   print "Energy = ", ham.H(0,0), " a.u."

    # Electronic DOFs
    el = Electronic(1,0)


    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)

    # Initialize MD variables
    nve_md.nve_md_init(syst, mol, el, ham)



    #=================== Propagation ====================
 
    ########################## Cooling #################################

    md = MD({"max_step":1,"ensemble":"NVT","integrator":"DLML","terec_exp_size":10,"dt":20.0,"n_medium":1,"n_fast":1,"n_outer":1})
    md.show_info()
    
    # Thermostat
    therm = Thermostat({"Temperature":278.0,"Q":100.0,"thermostat_type":"Nose-Hoover","nu_therm":0.001,"NHC_size":5})
    therm.show_info()


    ST = State() 
    ST.set_system(syst)
    ST.set_thermostat(therm)
    ST.set_md(md)
    ST.init_md(mol, el, ham, rnd)    



    f = open("_en_cooling.txt","w")
    f.close()


    md.dt = 20.0 
    md.max_step = 100
    for i in xrange(100):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_cooling.xyz",i)

        ST.run_md(mol, el, ham)
        ekin = ST.E_kin;    epot = ST.E_pot
        ST.cool()

        f = open("_en_cooling.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ekin, epot, ST.E_tot, ST.H_NP, ST.curr_T ))
        f.close()


    ########################## Production MD #################################

    syst.init_atom_velocities(300.0, rnd) 

    f = open("_en_md.txt","w")
    f.close()
    md.dt = 40.0 
    md.max_step = 10

    for i in xrange(1000):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)

        ST.run_md(mol, el, ham)

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ST.E_kin, ST.E_pot, ST.E_tot, ST.H_NP, ST.curr_T ))
        f.close()


main()

