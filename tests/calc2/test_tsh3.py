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

######################################################################
#
# Here we try to work on some acceleration of computations
#
######################################################################

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


ntraj = 1
use_boltz_factor = 0
T = 300.0
rnd = Random()
do_rescaling = 1 
do_reverse = 1
nu_therm = 0.01
rep = 1 # adiabatic
#ham_indx = 1
nel = 2
dt = 20.0
istate = 0

#=========== STEP 2:  Create system and load a molecule ================
SYST = []
for i in xrange(ntraj):
    syst = System()
#    Load_Molecule(U, syst, os.getcwd()+"/ch4.pdb", "pdb_1")
#    Load_Molecule(U, syst, os.getcwd()+"/c2h4.pdb", "pdb_1")
#    Load_Molecule(U, syst, os.getcwd()+"/bh.pdb", "pdb_1")
    LoadMolecule.Load_Molecule(U, syst, os.getcwd()+"/azobenz.pdb", "pdb_1")

    syst.init_fragments()
    syst.init_fragment_velocities(300.0 )
    # In this case Atoms are Fragments, so:
    print "syst.Number_of_atoms = ",syst.Number_of_atoms
    print "syst.Number_of_fragments = ",syst.Number_of_fragments


    for n in xrange(syst.Number_of_atoms):
        print n, syst.Atoms[n].Atom_RB.rb_p.x, syst.Atoms[n].Atom_RB.rb_p.y, syst.Atoms[n].Atom_RB.rb_p.z
        syst.Atoms[n].Atom_RB.rb_p = syst.Fragments[n].Group_RB.rb_p
        print n, syst.Atoms[n].Atom_RB.rb_p.x, syst.Atoms[n].Atom_RB.rb_p.y, syst.Atoms[n].Atom_RB.rb_p.z

    SYST.append(syst)



nnucl = 3*SYST[0].Number_of_atoms



#=========== STEP 3: Create control parameters (setting computation options) ================
# Creating Hamiltonian
HAM = []
for i in xrange(ntraj):
    ham = Hamiltonian_Atomistic(nel, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("QM")
    ham.set_system(SYST[i])
    ham.init_qm_Hamiltonian("control_parameters.dat")
    ham.set_rep(rep)


    # Create all excited states
    for e in range(1,nel):
        ham.add_excitation(0,1,e,1)
    HAM.append(ham)


therms = []  
for i in xrange(ntraj):
    therm = Thermostat({"nu_therm":nu_therm, "NHC_size":2, "Temperature":300.0, "thermostat_type":"Nose-Hoover" }) 
    therm.set_Nf_t(1);
    therm.set_Nf_r(0);
    therm.init_nhc();
    therms.append(therm)



MOL = []
for i in xrange(ntraj):
    mol = Nuclear(3*syst.Number_of_atoms)
    SYST[i].extract_atomic_q(mol.q)
    SYST[i].extract_atomic_p(mol.p)
    SYST[i].extract_atomic_f(mol.f)
    SYST[i].extract_atomic_mass(mol.mass)

    MOL.append(mol)


EL = []
for i in xrange(ntraj):
    el = Electronic(nel,istate)  # start on the 1-st excited state!!!
    EL.append(el)



t1 = Timer()
t2 = Timer()

f = open("relax.txt","w")
f.close()


# Initialization
Etot = []

for i in xrange(ntraj):
    v = []
    for n in xrange(nnucl):
        v.append(MOL[i].p[n]/MOL[i].mass[n])
    print v
    HAM[i].set_v( v )

    t2.start()
    epot = compute_forces(MOL[i], EL[i], HAM[i], 1) # 0 - MF, 1 - FSSH
    print "Time to compute epot = ", t2.stop()

    ekin = compute_kinetic_energy(MOL[i])
    Etot.append(epot+ekin)


# Propagation
t1.start()

print "Step 2: Propagation"
for snap in xrange(350):

    if snap==10:
        for i in xrange(ntraj):
            EL[i] = Electronic(nel,1)
            EL[i].istate = 1
              
    
    for t in xrange(1):


        for i in xrange(ntraj):

            # el-dyn
            EL[i].propagate_electronic(0.5*dt, HAM[i])



            #******************** Nuclear dynamics *****************************

            # Thermostat:            
            for n in xrange(nnucl):               
                MOL[i].p[n] = MOL[i].p[n] * therms[i].vel_scale(0.5*dt)


            MOL[i].propagate_p(0.5*dt)            
            MOL[i].propagate_q(dt)


            v = []
            for n in xrange(nnucl):
                v.append(MOL[i].p[n]/MOL[i].mass[n])
            HAM[i].set_v(v)


            t2.start()
            epot = compute_forces(MOL[i], EL[i], HAM[i], 1)
            ekin = compute_kinetic_energy(MOL[i])
            print "Time to compute forces = ", t2.stop()


            therms[i].propagate_nhc(dt, ekin, 0.0, 0.0)


            MOL[i].propagate_p(0.5*dt)


            # Thermostat
            for n in xrange(nnucl):               
                MOL[i].p[n] = MOL[i].p[n] * therms[i].vel_scale(0.5*dt)

            #******************** Nuclear dynamics ends **************************



            # el-dyn
            v = []
            for n in xrange(nnucl):
                v.append(MOL[i].p[n]/MOL[i].mass[n])
            HAM[i].set_v(v)
    
            EL[i].propagate_electronic(0.5*dt, HAM[i])

            ekin = compute_kinetic_energy(MOL[i])


            Etot[i] = epot + ekin        




            # Now, incorporate surface hop

            t2.start()

            g = MATRIX(nel, nel)
            compute_hopping_probabilities_fssh(MOL[i], EL[i], HAM[i], g, dt, use_boltz_factor, T) 
            print "Hopping probability matrix = "
            g.show_matrix()

            print "Time to compute hopping probabilities = ", t2.stop()

 
            t2.start()

            ksi = rnd.uniform(0.0, 1.0)
            EL[i].istate = hop(EL[i].istate, MOL[i], HAM[i], ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary

            print "Time to do actual hops = ", t2.stop()


    # Compute statistics for this snap    
    se_pops = [0.0] * nel
    sh_pops = [0.0] * nel

    for s in xrange(nel):
        for i in xrange(ntraj):
            se_pops[s] += EL[i].rho(s,s).real
        se_pops[s] = se_pops[s]/float(ntraj)

    for i in xrange(ntraj):
        sh_pops[ EL[i].istate ] += 1.0

    for s in xrange(nel):
        sh_pops[s] = sh_pops[s]/float(ntraj)

    E0  = HAM[0].Hvib(0,0).real
    E1  = HAM[0].Hvib(1,1).real
    E01 = HAM[0].Hvib(0,1).imag

    f = open("relax.txt","a")
    f.write("snap= %5i  se_pop(0)= %8.5f  se_pop(1)= %8.5f  sh_pop(0)= %8.5f sh_pop(1)= %8.5f Etot[0]= %8.5f E0= %8.5f E1= %8.5f E01= %8.5f\n" % (snap, se_pops[0], se_pops[1], sh_pops[0], sh_pops[1], Etot[0], E0, E1, E01 ))
    f.close()

    for i in xrange(ntraj):
        SYST[i].set_atomic_q(MOL[i].q)
        SYST[i].print_xyz("tsh_copy_"+str(i)+".xyz",snap)


print "time propagation is done"      
print "file is closed"
print "Time to complete = ", t1.stop(), " sec"


