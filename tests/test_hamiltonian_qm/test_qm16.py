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
# Tutorial: Atomistic Hamiltonian - excitonic excited state energies
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


#=========== STEP 1:  Create Universe and populate it ================
U = Universe()
LoadPT.Load_PT(U, "elements.dat", 1)


#=========== STEP 2:  Create system and load a molecule ================
syst = System()
#Load_Molecule(U, syst, os.getcwd()+"/c.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/c2.pdb", "pdb_1")
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/bh.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/co.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/ch4.pdb", "pdb_1")


print "Number of atoms in the system = ", syst.Number_of_atoms
atlst1 = range(0,syst.Number_of_atoms)


#=========== STEP 3: Create control parameters (setting computation options) ================


# Creating Hamiltonian
ham = Hamiltonian_Atomistic(3, 3*syst.Number_of_atoms)
ham.set_Hamiltonian_type("QM")
ham.set_system(syst)
ham.init_qm_Hamiltonian("control_parameters.dat")
ham.set_rep(1)

ham.add_excitation(0,1,1,1)
ham.add_excitation(0,1,2,1)

#ham.compute()
#print "Energy = ", ham.H(0,0), " a.u."
#print "Energy = ", ham.H(1,1), " a.u."
#print "Energy = ", ham.H(2,2), " a.u."

mol = Nuclear(3*syst.Number_of_atoms)
syst.extract_atomic_q(mol.q)


f = open("forces_ex.txt","w")
dx = 0.01
for n in xrange(50):
    mol.q[3] += 0.1  # x1
    ham.set_q(mol.q)
    ham.compute()
    E0, E1, E2 = ham.H(0,0).real, ham.H(1,1).real, ham.H(2,2).real

    for i in xrange(6):
        print -ham.dHdq(0,0,i).real, -ham.dHdq(1,1,i).real, -ham.dHdq(2,2,i).real

    f0, f1, f2 = -ham.dHdq(0,0,3).real, -ham.dHdq(1,1,3).real, -ham.dHdq(2,2,3).real


    mol.q[3] -= dx  # x1
    ham.set_q(mol.q)
    ham.compute()
    E0_l, E1_l, E2_l = ham.H(0,0).real, ham.H(1,1).real, ham.H(2,2).real


    mol.q[3] += 2.0*dx  # x1
    ham.set_q(mol.q)
    ham.compute()
    E0_r, E1_r, E2_r = ham.H(0,0).real, ham.H(1,1).real, ham.H(2,2).real

    mol.q[3] -= dx  # x1
    ham.set_q(mol.q)

    f0_n = -(E0_r-E0_l)/(2.0*dx)
    f1_n = -(E1_r-E1_l)/(2.0*dx)
    f2_n = -(E2_r-E2_l)/(2.0*dx)

    f.write("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n" % (mol.q[3],  E0, E1, E2,  f0, f1, f2,  f0_n, f1_n, f2_n) )
    
f.close()




#print "Force 1 = ", ham.dHdq(0,0,0), " a.u."
#print "Force 3 = ", ham.dHdq(0,0,3), " a.u."




