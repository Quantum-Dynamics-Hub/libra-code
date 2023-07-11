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
# Tutorial: Compute core INDO Hamiltonian, solve eigenvalue problem, 
#  get populations and density matrices
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
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/c.pdb", "pdb_1")

print "Number of atoms in the system = ", syst.Number_of_atoms
atlst1 = range(0,syst.Number_of_atoms)


#=========== STEP 3: Create control parameters (setting computation options) ================
prms = Control_Parameters()
#get_parameters_from_file("control_parameters_indo.dat", prms)
get_parameters_from_file("control_parameters_eht.dat", prms)
print "guess type = ", prms.guess_type


#=========== STEP 4:  Create model parameters and load them from file (using control parameters options) ================
modprms = Model_Parameters()

# Initialize/read model parameters (need basis info)
print "Setting parameters"
if(prms.hamiltonian=="eht"):
    set_parameters_eht(prms, modprms)
elif(prms.hamiltonian=="indo"):
    set_parameters_indo(prms, modprms);


#=========== STEP 5: Set basis (STO-3G_DZ) ================
Nelec = 0;
Norb = 0;

#------- Input --------------
mol_at_types = StringList()
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
x_period = 0
y_period = 0
z_period = 0
t1 = VECTOR()
t2 = VECTOR()
t3 = VECTOR()



update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, Sao)
print "AO overlap matrix"
Sao.show_matrix()

#=========== STEP 8: Parameters and core Hamiltonian ================
eri = doubleList()
V_AB = doubleList()
opt = 1  # 1 - for INDO, 0 - for CNDO/CNDO2
Hao = MATRIX(Norb, Norb)
debug = 0

     
if(prms.hamiltonian=="indo"):
    Sao.Init_Unit_Matrix(1.0)
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,0)
    Hamiltonian_core_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao,  Sao, debug)

elif(prms.hamiltonian=="eht"):
    Hamiltonian_core_eht(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao,  Sao, debug)

Hao.show_matrix()



#=========== STEP 9: Combine everything we need as an Electronic_Structure object ================
el = Electronic_Structure(Norb)
el.set_Hao(Hao)
el.set_Sao(Sao);

E_alp = MATRIX(Norb, Norb)
C_alp = MATRIX(Norb, Norb)
E_bet = MATRIX(Norb, Norb)
C_bet = MATRIX(Norb, Norb)


solve_eigen(Hao, Sao, E_alp, C_alp, 0)
solve_eigen(Hao, Sao, E_bet, C_bet, 0)




#
E_alp.show_matrix()

e_alp = []
e_bet = []
for i in xrange(Norb):
    e_alp.append(E_alp.get(i,i))
    e_bet.append(E_bet.get(i,i))

Nelec_alp = Nelec/2
Nelec_bet = Nelec - Nelec_alp
print "Nelec_alp = ", Nelec_alp
print "Nelec_bet = ", Nelec_bet


degen = 1.0
kT = 0.025
etol = 0.0001
Ef_alp = fermi_energy(e_alp, Nelec_alp, degen, kT, etol)
Ef_bet = fermi_energy(e_bet, Nelec_bet, degen, kT, etol)
print "Fermi energy (alpha) = ", Ef_alp
print "Fermi energy (beta) = ", Ef_bet

bnds_alp = order_bands(E_alp)
bnds_bet = order_bands(E_bet)
print "Printing orbital energies"
print bnds_alp
print bnds_bet

pop_opt = 0  #  0 -  integer populations,  1 - Fermi distribution              
occ_alp = populate_bands(Nelec_alp, degen, kT, etol, pop_opt, bnds_alp)
occ_bet = populate_bands(Nelec_bet, degen, kT, etol, pop_opt, bnds_bet)
print "Printing orbital occupations"
print occ_alp
print occ_bet


P_alp = compute_density_matrix(occ_alp, C_alp)
P_bet = compute_density_matrix(occ_bet, C_bet)

print "Density matrix (alpha)\n"
P_alp.show_matrix()
print "Density matrix (bet)\n"
P_bet.show_matrix()


