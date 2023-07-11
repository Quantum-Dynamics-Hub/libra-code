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
# Tutorial: Make a function to compute forces - tested - forces comply with F_ij = -F_ji
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



def force(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, Norb,c, x_period, y_period, z_period, t1, t2, t3):
    
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

    Hamiltonian_core_deriv_indo(syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, DF, c, dHao_dx, dHao_dy, dHao_dz, dSao_dx, dSao_dy, dSao_dz )
    Hamiltonian_Fock_derivs_indo(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, c, dHao_dx, dHao_dy, dHao_dz, dFao_alp_dx, dFao_alp_dy, dFao_alp_dz, dFao_bet_dx, dFao_bet_dy, dFao_bet_dz)
    update_derivative_coupling_matrix(x_period, y_period, z_period, t1, t2, t3, atom_to_ao_map, ao_to_atom_map, basis_ao, c, Dao_x, Dao_y, Dao_z);

    #print "Dao_x = ";    Dao_x.show_matrix()
    #print "Dao_y = ";    Dao_y.show_matrix()
    #print "Dao_z = ";    Dao_z.show_matrix()
    #print "Checking properties of Dao_x vs. dSao_dx"    
    #tmp = ( Dao_x + Dao_x.T() ) - dSao_dx
    #tmp.show_matrix()
    #print "Checking properties of Dao_y vs. dSao_dy"    
    #tmp = ( Dao_y + Dao_y.T() ) - dSao_dy
    #tmp.show_matrix()
    #print "Checking properties of Dao_z vs. dSao_dz"    
    #tmp = ( Dao_z + Dao_z.T() ) - dSao_dz
    #tmp.show_matrix()

    # Because in INDO  S = I => Dao = 0.0
    for i in xrange(Norb):
        for j in xrange(Norb):
            dSao_dx.set(i,j,0.0)
            dSao_dy.set(i,j,0.0)
            dSao_dz.set(i,j,0.0)
            Dao_x.set(i,j,0.0)
            Dao_y.set(i,j,0.0)
            Dao_z.set(i,j,0.0)


    C_alp = el.get_C_alp()
    C_bet = el.get_C_bet()
    E_alp = el.get_E_alp()
    E_bet = el.get_E_bet()

    T = E_alp * C_alp * Dao_x * C_alp.T()
    A_ax = C_alp * dFao_alp_dx * C_alp.T() - (T + T.T())
    T = E_alp * C_alp * Dao_y * C_alp.T()
    A_ay = C_alp * dFao_alp_dy * C_alp.T() - (T + T.T())
    T = E_alp * C_alp * Dao_z * C_alp.T()
    A_az = C_alp * dFao_alp_dz * C_alp.T() - (T + T.T())

    T = E_bet * C_bet * Dao_x * C_bet.T()
    A_bx = C_bet * dFao_bet_dx * C_bet.T() - (T + T.T())
    T = E_bet * C_bet * Dao_y * C_bet.T()
    A_by = C_bet * dFao_bet_dy * C_bet.T() - (T + T.T())
    T = E_bet * C_bet * Dao_z * C_bet.T()
    A_bz = C_bet * dFao_bet_dz * C_bet.T() - (T + T.T())

    dEa_dx = MATRIX(Norb,Norb)
    dEa_dy = MATRIX(Norb,Norb)
    dEa_dz = MATRIX(Norb,Norb)
    dEb_dx = MATRIX(Norb,Norb)
    dEb_dy = MATRIX(Norb,Norb)
    dEb_dz = MATRIX(Norb,Norb)

    Dmo_a_x = MATRIX(Norb,Norb)
    Dmo_a_y = MATRIX(Norb,Norb)
    Dmo_a_z = MATRIX(Norb,Norb)
    Dmo_b_x = MATRIX(Norb,Norb)
    Dmo_b_y = MATRIX(Norb,Norb)
    Dmo_b_z = MATRIX(Norb,Norb)


    for i in xrange(Norb):
        dEa_dx.set(i,i,A_ax.get(i,i))
        dEa_dy.set(i,i,A_ay.get(i,i))
        dEa_dz.set(i,i,A_az.get(i,i))
        dEb_dx.set(i,i,A_bx.get(i,i))
        dEb_dy.set(i,i,A_by.get(i,i))
        dEb_dz.set(i,i,A_bz.get(i,i))

        for j in xrange(Norb):
            if(i!=j):
                Dmo_a_x.set(i,j,A_ax.get(i,j)/(E_alp.get(j,j)-E_alp.get(i,i)))
                Dmo_a_y.set(i,j,A_ay.get(i,j)/(E_alp.get(j,j)-E_alp.get(i,i)))
                Dmo_a_z.set(i,j,A_az.get(i,j)/(E_alp.get(j,j)-E_alp.get(i,i)))
                Dmo_b_x.set(i,j,A_bx.get(i,j)/(E_bet.get(j,j)-E_bet.get(i,i)))
                Dmo_b_y.set(i,j,A_by.get(i,j)/(E_bet.get(j,j)-E_bet.get(i,i)))
                Dmo_b_z.set(i,j,A_bz.get(i,j)/(E_bet.get(j,j)-E_bet.get(i,i)))

#    print "dEa_dx = "
#    dEa_dx.show_matrix()
#    print "dEa_dy = "
#    dEa_dy.show_matrix()
#    print "dEa_dz = "
#    dEa_dz.show_matrix()
#    print "dEb_dx = "
#    dEb_dx.show_matrix()
#    print "dEb_dy = "
#    dEb_dy.show_matrix()
#    print "dEb_dz = "
#    dEb_dz.show_matrix()

#    print "Dmo_a_x = "
#    Dmo_a_x.show_matrix()
#    print "Dmo_a_y = "
#    Dmo_a_y.show_matrix()
#    print "Dmo_a_z = "
#    Dmo_a_z.show_matrix()
#    print "Dmo_b_x = "
#    Dmo_b_x.show_matrix()
#    print "Dmo_b_y = "
#    Dmo_b_y.show_matrix()
#    print "Dmo_b_z = "
#    Dmo_b_z.show_matrix()


    # P = C * O * C.T(), and C^T * S * C = O  =>  O = C^T * S * P * S * C
    Oa = C_alp.T() * Sao * el.get_P_alp() * Sao * C_alp
    Ob = C_bet.T() * Sao * el.get_P_bet() * Sao * C_bet

    T = C_alp * C_alp.T() * Dao_x * C_alp
    dP_alp_dx = C_alp * (Dmo_a_x * Oa - Oa * Dmo_a_x) * C_alp.T() - (T + T.T())

    T = C_alp * C_alp.T() * Dao_y * C_alp
    dP_alp_dy = C_alp * (Dmo_a_y * Oa - Oa * Dmo_a_y) * C_alp.T() - (T + T.T())

    T = C_alp * C_alp.T() * Dao_z * C_alp
    dP_alp_dz = C_alp * (Dmo_a_z * Oa - Oa * Dmo_a_z) * C_alp.T() - (T + T.T())

    T = C_bet * C_bet.T() * Dao_x * C_bet
    dP_bet_dx = C_bet * (Dmo_b_x * Ob - Ob * Dmo_b_x) * C_bet.T() - (T + T.T())

    T = C_bet * C_bet.T() * Dao_y * C_bet
    dP_bet_dy = C_bet * (Dmo_b_y * Ob - Ob * Dmo_b_y) * C_bet.T() - (T + T.T())

    T = C_bet * C_bet.T() * Dao_z * C_bet
    dP_bet_dz = C_bet * (Dmo_b_z * Ob - Ob * Dmo_b_z) * C_bet.T() - (T + T.T())


    F_x =       ( dP_alp_dx * (Hao + el.get_Fao_alp()) + el.get_P_alp() * (dHao_dx + dFao_alp_dx) ).tr()
    F_x = F_x + ( dP_bet_dx * (Hao + el.get_Fao_bet()) + el.get_P_bet() * (dHao_dx + dFao_bet_dx) ).tr()

    F_y =       ( dP_alp_dy * (Hao + el.get_Fao_alp()) + el.get_P_alp() * (dHao_dy + dFao_alp_dy) ).tr()
    F_y = F_y + ( dP_bet_dy * (Hao + el.get_Fao_bet()) + el.get_P_bet() * (dHao_dy + dFao_bet_dy) ).tr()

    F_z =       ( dP_alp_dz * (Hao + el.get_Fao_alp()) + el.get_P_alp() * (dHao_dz + dFao_alp_dz) ).tr()
    F_z = F_z + ( dP_bet_dz * (Hao + el.get_Fao_bet()) + el.get_P_bet() * (dHao_dz + dFao_bet_dz) ).tr()

    print " Force = ", F_x, F_y, F_x

    return [F_x, F_y, F_x]




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

if(prms.hamiltonian=="eht" or prms.hamiltonian=="geht" or prms.hamiltonian=="geht1" or prms.hamiltonian=="geht2"):
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


print "Orbitals"
for i in xrange(Norb):
    print "Orbital", i, basis_ao[i].show_info()

c = 0
F = force(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, Norb, c, x_period, y_period, z_period, t1, t2, t3)
print F

c = 1
F = force(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map, Hao, Sao, Norb, c, x_period, y_period, z_period, t1, t2, t3)
print F









