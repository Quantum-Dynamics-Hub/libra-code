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
import math
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


#=========== STEP 1:  Create Universe and populate it ================
U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.dat", 1)


#=========== STEP 2:  Create system and load a molecule ================
syst = System()
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/azulene.pdb", "iqmol_pdb")
#Load_Molecule(U, syst, os.getcwd()+"/benzene.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/c.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/c2.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/ch4.pdb", "pdb_1")

print "Number of atoms in the system = ", syst.Number_of_atoms
atlst1 = range(0,syst.Number_of_atoms)


#=========== STEP 3: Create control parameters (setting computation options) ================
prms = Control_Parameters()
#get_parameters_from_file("control_parameters_indo.dat", prms)
get_parameters_from_file("control_parameters_eht.dat", prms)

#=========== STEP 4:  Create model parameters and load them from file (using control parameters options) ================
modprms = Model_Parameters()

# Initialize/read model parameters (need basis info)
print "Setting parameters"
if(prms.hamiltonian=="eht"):
    set_parameters_eht(prms, modprms)
    modprms.eht_k.show()
elif(prms.hamiltonian=="indo"):
    set_parameters_indo(prms, modprms)


#=========== STEP 5: Set basis (STO-3G_DZ) ================
Nelec, Norb = 0, 0

#------- Input --------------
mol_at_types = StringList();
R = VECTORList()

for i in xrange(syst.Number_of_atoms):
    mol_at_types.append(syst.Atoms[i].Atom_element)
    R.append(syst.Atoms[i].Atom_RB.rb_cm)

#-------- Output -----------
basis = AOList()
atom_to_ao_map = intMap()
ao_to_atom_map = intList()

verb = 0
basis_ao, Nelec, Norb, atom_to_ao_map, ao_to_atom_map = set_basis_STO_3G_DZ(mol_at_types, R, modprms, verb)


#=========== STEP 6: Depending on hamiltonian to use, set internal parameters ================

if(prms.hamiltonian=="eht"):
    set_parameters_eht_mapping(modprms, basis_ao)    
    set_parameters_eht_mapping1(modprms,syst.Number_of_atoms,mol_at_types)


#=========== STEP 7: Overlap matrix ================

Sao = MATRIX(Norb, Norb)
x_period = 0; y_period = 0; z_period = 0
t1 = VECTOR(); t2 = VECTOR(); t3 = VECTOR()

update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, Sao);
print "AO overlap matrix"
Sao.show_matrix()

#=========== STEP 8: Parameters ================
eri = doubleList()
V_AB = doubleList()
opt = 1  # 1 - for INDO, 0 - for CNDO/CNDO2
Hao = MATRIX(Norb, Norb)
debug = 1
     
if(prms.hamiltonian=="indo"):
    Sao.Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);
    Hamiltonian_core(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao,  Sao, debug)

elif(prms.hamiltonian=="eht"):
    Hamiltonian_core_eht(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao,  Sao, debug)

print "Core Hamiltonian"
Hao.show_matrix()


Nelec_alp = Nelec/2
Nelec_bet = Nelec - Nelec_alp
print "Nelec_alp = ", Nelec_alp
print "Nelec_bet = ", Nelec_bet

degen = 1.0
kT = 0.025
etol = 0.0001
pop_opt = 0  #  0 -  integer populations,  1 - Fermi distribution              

res_alp = Fock_to_P(Hao, Sao, Nelec_alp, degen, kT, etol, pop_opt)
res_bet = Fock_to_P(Hao, Sao, Nelec_bet, degen, kT, etol, pop_opt)


print "Eigenvalues (alp):\n"
res_alp[0].show_matrix()
print "Eigenvalues (bet):\n"
res_bet[0].show_matrix()

print "Eigenvectors(alp):\n"
res_alp[1].show_matrix()
print "Eigenvectors(bet):\n"
res_bet[1].show_matrix()


print "Density matrix(alp):\n"
res_alp[2].show_matrix()
print "Density matrix(bet):\n"
res_bet[2].show_matrix()


print "Bands(alp):\n"
print res_alp[3]
print "Bands(bet):\n"
print res_bet[3]


print "Occupations(alp):\n"
print res_alp[4]
print "Occupations(bet):\n"
print res_bet[4]



el = Electronic_Structure(Norb)
el.Nocc_alp = Nelec_alp
el.Nocc_bet = Nelec_bet
el.set_Hao(Hao)
el.set_Sao(Sao)
el.set_P_alp(res_alp[2])
el.set_P_bet(res_bet[2])


if(prms.hamiltonian=="indo"):
    Hamiltonian_Fock_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map)
elif(prms.hamiltonian=="eht"):
    Hamiltonian_Fock_eht(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map)


print "Fock matrix at first iteration (alp)"
el.get_Fao_alp().show_matrix()

#===============  Now to SCF iterations =======================

E = scf(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, 0); 


degen = 1.0
kT = 0.025
etol = 0.0001
pop_opt = 0  #  0 -  integer populations,  1 - Fermi distribution              

res_alp = Fock_to_P(el.get_Fao_alp(), el.get_Sao(), Nelec_alp, degen, kT, etol, pop_opt)
res_bet = Fock_to_P(el.get_Fao_bet(), el.get_Sao(), Nelec_bet, degen, kT, etol, pop_opt)


print "Bands(alp)    Occupations(alp)       Bands(bet)    Occupations(bet)"
for j in xrange(Norb):
     print "%5.3i %12.8f   %12.8f  %12.8f   %12.8f" % (j, res_alp[3][j][1], res_alp[4][j][1], res_bet[3][j][1], res_bet[4][j][1])


vlist = intList()    
vlist.append(23)
vlist.append(24)

prms.orbs = vlist
prms.nx_grid = 40
prms.ny_grid = 40
prms.nz_grid = 40
prms.charge_density_prefix = "char_dens/"

charge_density( el, syst, basis_ao, prms)




