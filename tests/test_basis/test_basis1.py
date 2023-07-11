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
# Tutorial: Here we manually construct AO (STO-3G_DZ) basis for H2O molecule
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
LoadPT.Load_PT(U, "elements.dat", 0)


#=========== STEP 2:  Create system and load a molecule ================
syst = System()
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/1water.ent", "pdb")

print "Number of atoms in the system = ", syst.Number_of_atoms
atlst1 = range(0,syst.Number_of_atoms)


#=========== STEP 3: Create control parameters (setting computation options) ================
prms = Control_Parameters()
get_parameters_from_file("control_parameters.dat", prms)
print "guess type = ", prms.guess_type


#=========== STEP 4:  Create model parameters and load them from file (using control parameters options) ================
modprms = Model_Parameters()

# Initialize/read model parameters (need basis info)
print "Setting parameters"
if(prms.hamiltonian=="eht" or prms.hamiltonian=="geht"):
    set_parameters_eht(prms, modprms)
elif(prms.hamiltonian=="indo"):
    set_parameters_indo(prms, modprms);
elif(prms.hamiltonian=="geht1"):
    set_parameters_geht1(prms, modprms); 
elif(prms.hamiltonian=="geht2"):
    set_parameters_geht2(prms, modprms); 


#=========== STEP 5: Set basis (STO-3G_DZ) ================
Nelec = 0;
Norb = 0;

print "Setting STO-3G_DZ basis:";
basis_ao = []

for i in xrange(syst.Number_of_atoms):
    print syst.Atoms[i].Atom_element, len(modprms.PT[syst.Atoms[i].Atom_element].IP)

    Nelec = Nelec + modprms.PT[syst.Atoms[i].Atom_element].Nval
    print "Atom", syst.Atoms[i].Atom_element, "contributes", modprms.PT[syst.Atoms[i].Atom_element].Nval, "electrons"

    for shell in modprms.PT[syst.Atoms[i].Atom_element].IP.keys():
        print shell, modprms.PT[syst.Atoms[i].Atom_element].IP[shell]

        Nzeta = modprms.PT[syst.Atoms[i].Atom_element].Nzeta[shell]
        Nquant = modprms.PT[syst.Atoms[i].Atom_element].Nquant[shell]
        IP = modprms.PT[syst.Atoms[i].Atom_element].IP[shell]
        exp1 = modprms.PT[syst.Atoms[i].Atom_element].zetas[shell][0]
        exp2 = modprms.PT[syst.Atoms[i].Atom_element].zetas[shell][1]
        coeff1 = modprms.PT[syst.Atoms[i].Atom_element].coeffs[shell][0]
        coeff2 = modprms.PT[syst.Atoms[i].Atom_element].coeffs[shell][1]


        add_basis_ao(syst.Atoms[i].Atom_element, syst.Atoms[i].Atom_RB.rb_cm,
                   shell[1], Nzeta, Nquant, IP, exp1, exp2, coeff1, coeff2, basis_ao)

    
        print len(basis_ao), "orbitals added"
        print Nzeta, Nquant, IP, exp1, exp2, coeff1, coeff2


Norb = len(basis_ao)
print "Total number of orbitals added is = ", Norb
print "Total number of electrons is = ", Nelec

# Print info about loaded AOs
for i in xrange(Norb):
    print "orbital ",i
    basis_ao[i].show_info()


