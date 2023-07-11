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
# Compute derivatives of Sao, Hao, Fao, and also Dao in AO basis
# Now, get step-by-step to the dE/dR, Dmo and forces (only for alpha component)
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
prms = Control_Parameters()
get_parameters_from_file("control_parameters_indo.dat", prms)
#get_parameters_from_file("control_parameters_eht.dat", prms)
print "guess type = ", prms.guess_type


#=========== STEP 4:  Create model parameters and load them from file (using control parameters options) ================
modprms = Model_Parameters()

# Initialize/read model parameters (need basis info)
print "Setting parameters"
if(prms.hamiltonian=="eht"):
    set_parameters_eht(prms, modprms)
elif(prms.hamiltonian=="indo"):
    set_parameters_indo(prms, modprms)



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



update_overlap_matrix(x_period, y_period, z_period, t1, t2, t3, basis_ao, Sao);
print "AO overlap matrix"
Sao.show_matrix()

#=========== STEP 8: Parameters ================
eri = doubleList()
V_AB = doubleList()
opt = 1  # 1 - for INDO, 0 - for CNDO/CNDO2


if(prms.hamiltonian=="indo"):
    Sao.Init_Unit_Matrix(1.0);  
    indo_core_parameters(syst, basis_ao, modprms, atom_to_ao_map, ao_to_atom_map, opt,1);



Hao = MATRIX(Norb, Norb)
debug = 1
Hamiltonian_core(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao,  Sao, debug)
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


Hamiltonian_Fock(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map)

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
     print "%12.8f   %12.8f  %12.8f   %12.8f" %(res_alp[3][j][1], res_alp[4][j][1], res_bet[3][j][1], res_bet[4][j][1])



    
dHao_dx = MATRIX(Norb, Norb)
dHao_dy = MATRIX(Norb, Norb)
dHao_dz = MATRIX(Norb, Norb)
dSao_dx = MATRIX(Norb, Norb)
dSao_dy = MATRIX(Norb, Norb)
dSao_dz = MATRIX(Norb, Norb)
dFao_alp_dx = MATRIX(Norb, Norb)
dFao_alp_dy = MATRIX(Norb, Norb)
dFao_alp_dz = MATRIX(Norb, Norb)
dFao_bet_dx = MATRIX(Norb, Norb)
dFao_bet_dy = MATRIX(Norb, Norb)
dFao_bet_dz = MATRIX(Norb, Norb)
Dao_x = MATRIX(Norb, Norb)
Dao_y = MATRIX(Norb, Norb)
Dao_z = MATRIX(Norb, Norb)


DF = 0

for c in [0,1]:
    print "Starting core derivatives"

    if(prms.hamiltonian=="indo"):
        Hamiltonian_core_deriv_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF, c, dHao_dx, dHao_dy, dHao_dz, dSao_dx, dSao_dy, dSao_dz )
        Hamiltonian_Fock_derivs_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, c, dHao_dx, dHao_dy, dHao_dz, dFao_alp_dx, dFao_alp_dy, dFao_alp_dz, dFao_bet_dx, dFao_bet_dy, dFao_bet_dz)
    elif (prms.hamiltonian=="eht"):
        print "Core matrix derivatives are not yet implemented for EHT..."
        sys.exit(0)

    update_derivative_coupling_matrix(x_period, y_period, z_period, t1, t2, t3, atom_to_ao_map, ao_to_atom_map, basis_ao, c, Dao_x, Dao_y, Dao_z);


    # Here goes the funny stuff
    if(prms.hamiltonian=="indo"):
        # Use a simplified version (orbitals are assumed to be orthonormal, the overlap is identity)

    else:

        C_alp = MATRIX(Norb, Norb);  C_alp = el.get_C_alp()
        E_alp = el.get_E_alp()

        T = MATRIX(Norb, Norb);    T = E_alp * C_alp * Dao_x * C_alp.T()

        A = MATRIX(Norb,Norb)
        A = C_alp * dFao_alp_dx * C_alp.T() - (T + T.T())
        A.show_matrix()

        dE_dx = MATRIX(Norb,Norb)
        Dmo_x = MATRIX(Norb,Norb)

    for i in xrange(Norb):
        dE_dx.set(i,i,A.get(i,i))
        for j in xrange(Norb):
            if(i!=j):
                Dmo_x.set(i,j,A.get(i,j)/(E_alp.get(j,j)-E_alp.get(i,i)))

    print "dE_dx["+str(c)+"]";   dE_dx.show_matrix()
    print "Dmo_x["+str(c)+"]";   Dmo_x.show_matrix()


    O = MATRIX(Norb,Norb)
    # P = C * O * C.T(), and C^T * S * C = I  =>  O = C^T * S * P * S * C
    O = C_alp.T() * Sao * el.get_P_alp() * Sao * C_alp
    O.show_matrix()

    T = C_alp * C_alp.T() * Dao_x * C_alp
    dP_alp_dx = C_alp * (Dmo_x * O - O * Dmo_x) * C_alp.T() - (T + T.T())

    print "dP_alp_dx["+str(c)+"] =";  dP_alp_dx.show_matrix()

    F_x = ( el.get_P_alp() * (dHao_dx + dFao_alp_dx) ).tr()   
    F_x = F_x + ( el.get_P_bet() * (dHao_dx + dFao_bet_dx) ).tr() 
    F_x = -0.5*F_x

    print " Force["+str(c)+"] = ", F_x


    # checking:
    print "Checking properties of Dao vs. dSao"
    tmp = MATRIX(Norb, Norb)
    tmp = ( Dao_x + Dao_x.T() ) - dSao_dx
    tmp.show_matrix()







