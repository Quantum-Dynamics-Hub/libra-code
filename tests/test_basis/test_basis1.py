#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
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
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/chemobjects")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_QM/Model_Parameters")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/qchem/qobjects")
sys.path.insert(1,cwd+"/../../_build/src/qchem/basis")
sys.path.insert(1,cwd+"/../../_build/src/converters")

print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygchemobjects import *
from cyghamiltonian import *
#from cyghamiltonian_atomistic import *
#from cyghamiltonian_qm import *
from cygcontrol_parameters import *
from cygmodel_parameters import *

from cygdyn import *

from cygqobjects import *
from cygbasis import *

from cygconverters import *


from LoadPT import * # Load_PT
from LoadMolecule import * # Load_Molecule



# Create Universe and populate it
U = Universe()
Load_PT(U, "elements.dat", 0)

syst = System()
Load_Molecule(U, syst, os.getcwd()+"/1water.ent", "pdb")


print "Number of atoms in the system = ", syst.Number_of_atoms
atlst1 = range(0,syst.Number_of_atoms)






prms = Control_Parameters()
get_parameters_from_file("control_parameters.dat", prms)

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



basis_ao = AOList()


# Set basis:
Nelec = 0;
Norb = 0;

print "Setting STO-3G_DZ basis:";
for i in xrange(syst.Number_of_atoms):
    print syst.Atoms[i].Atom_element, len(modprms.PT[syst.Atoms[i].Atom_element].IP)

#    print modprms.PT["O"]
    for shell in modprms.PT[syst.Atoms[i].Atom_element].IP.keys():
        print shell, modprms.PT[syst.Atoms[i].Atom_element].IP[shell]

        Nzeta = modprms.PT[syst.Atoms[i].Atom_element].Nzeta[shell]
        Nquant = modprms.PT[syst.Atoms[i].Atom_element].Nquant[shell]
        IP = modprms.PT[syst.Atoms[i].Atom_element].IP[shell]
        exp1 = modprms.PT[syst.Atoms[i].Atom_element].zetas[shell][0]
        exp2 = modprms.PT[syst.Atoms[i].Atom_element].zetas[shell][1]
        coeff1 = modprms.PT[syst.Atoms[i].Atom_element].coeffs[shell][0]
        coeff2 = modprms.PT[syst.Atoms[i].Atom_element].coeffs[shell][1]

        orbs = []

        add_basis_ao(syst.Atoms[i].Atom_element, syst.Atoms[i].Atom_RB.rb_cm,
                   shell, Nzeta, Nquant, IP, exp1, exp2, coeff1, coeff2, orbs)


        print len(orbs), "orbitals added"
        print Nzeta, Nquant, IP, exp1, exp2, coeff1, coeff2

#    // Compute the number of valence electrons
#    //Nelec += num_valence_elec(mol.Z[i]);
#    Nelec += modpar.PT[mol.at_type[i]].Nval;
#    cout<<"Atom "<<i<<" contributes "<<modpar.PT[mol.at_type[i]].Nval<<" electrons\n";

#  }// for i - all atoms in the system

#//  Now print all AOs:
#  cout<<"Printing atomic basis: \n";
#  for(i=0;i<Norb;i++){
#    cout<<i<<"  ao_name= "<<basis_ao[i].ao_name<<"  at_indx= "<<basis_ao[i].at_indx<<"  element= "<<basis_ao[i].element<<endl;
#  }







#  // Initialize nuclear subsystem
#  for(n=0;n<atoms.size();n++){
#
#    std::string elt = atoms[n].atom_type;
#    double zeff = modprms.PT[elt].Zeff;
#
#    mol.add_atom(elt, elt, modprms.PT[elt].Z, zeff, modprms.PT[elt].mass, atoms[n].R);
#
#  }// for n



#  //----------------------------------------------------------------------------------
#  //------------ Initialize basis of AOs and set atom-orbital mapping ----------------
#  set_basis_STO_3G_DZ(mol, modprms, basis_ao, Nelec, Norb);


#  if(prms.hamiltonian=="eht"||prms.hamiltonian=="geht"||prms.hamiltonian=="geht1"||prms.hamiltonian=="geht2"){  
#    set_parameters_eht_mapping(modprms,basis_ao); 
#    set_parameters_eht_mapping1(modprms,mol); 

#/*
#    cout<<"K matrix:\n";
#    for(i=0;i<basis_ao.size();i++){
#      for(j=0;j<basis_ao.size();j++){
#        cout<<modprms.meht_k.get_K_value(i,j)<<" ";
#      }
#      cout<<"\n";
#    }
#*/


print "guess type = ", prms.guess_type




#for i in atlist1:



# Creating Hamiltonian
#    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
#    ham.set_Hamiltonian_type("QM")

   
#    ham.set_system(syst)
#    ham.compute()


#    print "Energy = ", ham.H(0,0), " a.u."
#    print "Force 1 = ", ham.dHdq(0,0,0), " a.u."
#    print "Force 3 = ", ham.dHdq(0,0,3), " a.u."




