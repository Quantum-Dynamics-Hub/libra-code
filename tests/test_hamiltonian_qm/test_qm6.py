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
# Tutorial: SCF computations are hidden - use built-in function
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
ham = Hamiltonian_Atomistic(2, 3*syst.Number_of_atoms)
ham.set_Hamiltonian_type("QM")
ham.set_system(syst)
ham.init_qm_Hamiltonian("control_parameters_indo.dat")
#ham.init_qm_Hamiltonian("control_parameters_eht.dat")
ham.set_rep(1)
ham.add_excitation(0,1,1,1)


ham.compute()



for i in xrange(2):
    print "Energy = ", ham.H(i,i), " a.u."
print "Force 1 = ", -ham.dHdq(0,0,0), " a.u."
print "Force 3 = ", -ham.dHdq(0,0,3), " a.u."





#degen = 1.0
#kT = 0.025
#etol = 0.0001
#pop_opt = 0  #  0 -  integer populations,  1 - Fermi distribution              

#res_alp = Fock_to_P(el.get_Fao_alp(), el.get_Sao(), Nelec_alp, degen, kT, etol, pop_opt)
#res_bet = Fock_to_P(el.get_Fao_alp(), el.get_Sao(), Nelec_bet, degen, kT, etol, pop_opt)


#print "Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
#for j in xrange(Norb):
#     print "%12.8f   %12.8f  %12.8f   %12.8f" %(res_alp[3][j][1], res_alp[4][j][1], res_bet[3][j][1], res_bet[4][j][1])

    






