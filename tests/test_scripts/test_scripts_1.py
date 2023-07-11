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
# Tutorial: All Python-written procedures can be hidden in C++ routines
# Here we demonstrate md_init and md_run functions as the counterparts
# to Python-written scripts (see test_hamiltonian_mm/test_mm6*.py)
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

    print "Number of atoms in the system = ", syst.Number_of_atoms
    atlst1 = range(1,syst.Number_of_atoms+1)


    # Creating Hamiltonian
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)  
    ham.show_interactions_statistics()   
    ham.set_system(syst)
    ham.compute()
    print "Energy = ", ham.H(0,0), " a.u."
    print "Force 1 = ", ham.dHdq(0,0,0), " a.u."
    print "Force 3 = ", ham.dHdq(0,0,3), " a.u."


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


    md = MD({"integrator":"DLML", "dt":1.0,"ensemble":"NVE","max_step":10})
    therm = Thermostat({"thermostat_type":"Nose-Hoover", "NHC_size":2, "nu_therm":0.001, "Q":100.0, "Temperature":300.0})

    st = State()
    st.set_system(syst)
    #st.set_thermostat(therm)
    st.set_md(md)



    st.init_md(mol, el, ham)

    ########################## Cooling #################################

    f = open("_en_cooling.txt","w")
    f.close()


    for i in xrange(100):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_cooling.xyz",i)

        st.run_md(mol, el, ham)
        syst.cool()

        f = open("_en_cooling.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, st.E_kin, st.E_pot, st.E_tot))
        f.close()


    ########################## Production MD #################################

    f = open("_en_md.txt","w")
    f.close()
    md.dt = 20.0 

    for i in xrange(100):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)

        st.run_md(mol, el, ham)

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, st.E_kin, st.E_pot, st.E_tot))
        f.close()


