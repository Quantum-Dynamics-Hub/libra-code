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
# Tutorial: Now we are ready for adiabatic MD
###################################################################

import os
import sys
import math
import copy

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
#Load_Molecule(U, syst, os.getcwd()+"/bh.pdb", "pdb_1")
#Load_Molecule(U, syst, os.getcwd()+"/co.pdb", "pdb_1")
LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/ch4.pdb", "pdb_1")




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

el = Electronic_Structure(Norb)

el.Nocc_alp = Nelec/2          
el.Nocc_bet = Nelec - el.Nocc_alp
print "Nelec= ",Nelec, " Nocc_alp= ",el.Nocc_alp, "Nocc_bet= ",el.Nocc_bet

Epot = energy_and_forces(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map)

#sys.exit(0)

dt = 20.0  # 
t = 0.0

f = open("bh_opt.xyz","w")
f.close()

f = open("opt.txt","w")
f.close()

for n in xrange(syst.Number_of_atoms):
    syst.Atoms[n].Atom_RB.rb_p = VECTOR(0.0,0.0,0.0)
    print n, syst.Atoms[n].Atom_RB.rb_mass


syst.init_fragments()
syst.show_fragments()

####################### Cooling down ###################################
for step in xrange(50):

    for n in xrange(syst.Number_of_atoms):
        syst.Atoms[n].Atom_RB.rb_p += 0.5*dt * syst.Atoms[n].Atom_RB.rb_force
        syst.Atoms[n].Atom_RB.rb_cm += (dt/syst.Atoms[n].Atom_RB.rb_mass) * syst.Atoms[n].Atom_RB.rb_p

    Epot = energy_and_forces(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map)
    Ekin = 0.0
    for n in xrange(syst.Number_of_atoms):
        syst.Atoms[n].Atom_RB.rb_p += 0.5*dt * syst.Atoms[n].Atom_RB.rb_force 
        Ekin = Ekin + syst.Atoms[n].Atom_RB.ekin_tr()

    Etot = Ekin + Epot
    t = t + dt

    if step % 25 == 0:
        for n in xrange(syst.Number_of_atoms):
            syst.Atoms[n].Atom_RB.rb_p.x = 0.0
            syst.Atoms[n].Atom_RB.rb_p.y = 0.0
            syst.Atoms[n].Atom_RB.rb_p.z = 0.0


    f = open("opt.txt","a")
    f.write("t= %12.8f Ekin= %12.8f Epot= %12.8f Etot= %12.8f \n" % (t, Ekin, Epot, Etot) )
    f.close()
    syst.print_xyz("bh_opt.xyz",step)


####################### MD ###################################

#sys.exit(0)
dt = 20.0

f = open("bh.xyz","w")
f.close()

f = open("md.txt","w")
f.close()

t = 0.0

for step in xrange(250):
    for n in xrange(syst.Number_of_atoms):
        syst.Atoms[n].Atom_RB.rb_p = syst.Atoms[n].Atom_RB.rb_p + 0.5*dt * syst.Atoms[n].Atom_RB.rb_force
        syst.Atoms[n].Atom_RB.rb_cm = syst.Atoms[n].Atom_RB.rb_cm + dt* syst.Atoms[n].Atom_RB.rb_p/syst.Atoms[n].Atom_RB.rb_mass

#    Epot = energy_and_forces(syst, basis_ao, Nelec, Norb, atom_to_ao_map, ao_to_atom_map)
    Epot = energy_and_forces(el, syst, basis_ao, prms, modprms, atom_to_ao_map, ao_to_atom_map)

    Ekin = 0.0
    for n in xrange(syst.Number_of_atoms):
        syst.Atoms[n].Atom_RB.rb_p = syst.Atoms[n].Atom_RB.rb_p + 0.5*dt * syst.Atoms[n].Atom_RB.rb_force 
        Ekin = Ekin + syst.Atoms[n].Atom_RB.ekin_tr()

    #Epot = res[0]
    Etot = Ekin + Epot
    t = t + dt

    f = open("md.txt","a")
    f.write("t= %12.8f Ekin= %12.8f Epot= %12.8f Etot= %12.8f \n" % (t, Ekin, Epot, Etot) )
    f.close()
    
    syst.print_xyz("bh.xyz",step)


    





