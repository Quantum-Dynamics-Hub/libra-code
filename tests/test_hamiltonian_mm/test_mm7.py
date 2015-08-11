###################################################################
# Tutorial: Here is another version of MD: RB-MD via Systems objects
# also add thermostat & barostat - at some point, I'll have to complete this
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
sys.path.insert(1,cwd+"/../../_build/src/dyn")
#sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/hamiltonian_atomistic")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygchemobjects import *
from cyghamiltonian import *
from cygdyn import *
#from cyghamiltonian_atomistic import *


from LoadPT import * # Load_PT
from LoadMolecule import * # Load_Molecule
from LoadUFF import*



##############################################################

def init_md(syst, mol, el, ham):
    syst.set_atomic_q(mol.q)  # mol -> syst

    E_kin = 0.0
    for n in xrange(syst.Number_of_fragments):
        E_kin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    syst.zero_forces_and_torques()

    # Compute forces and 
    syst.extract_atomic_q(mol.q) # syst -> mol
    E_pot = compute_potential_energy(mol, el, ham, 1)  # 1 - FSSH forces
    compute_forces(mol, el, ham, 1)
    syst.set_atomic_f(mol.f)  # mol -> syst
    syst.update_fragment_forces_and_torques();

    E_tot = E_pot + E_kin;

    return [E_kin, E_pot, E_tot]


def md_step(syst, mol, el, ham, dt, integrator):

    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    ps = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Propagate translational DOFs 
        syst.Fragments[n].Group_RB.shift_position(dt * syst.Fragments[n].Group_RB.rb_p * syst.Fragments[n].Group_RB.rb_iM);
 
        # Propagate rotational DOFs
        if integrator=="Jacobi":
            syst.Fragments[n].Group_RB.propagate_exact_rb(dt)
        elif integrator=="DLML":
            ps = syst.Fragments[n].Group_RB.propagate_dlml(dt)

    # Update atomic coordinates
    for n in xrange(syst.Number_of_fragments):
        syst.update_atoms_for_fragment(n)

    #=========== Potential ===================
    syst.zero_forces_and_torques()
    syst.extract_atomic_q(mol.q)

    epot = compute_potential_energy(mol, el, ham, 1)  # 1 - FSSH forces
    compute_forces(mol, el, ham, 1)

    syst.set_atomic_f(mol.f)
    syst.update_fragment_forces_and_torques()


    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    etot = ekin+epot
    return [ekin, epot, etot]


def ekin_tr(syst):
    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        ekin += syst.Fragments[n].Group_RB.ekin_tr() 
    return ekin

def ekin_rot(syst):
    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        ekin += syst.Fragments[n].Group_RB.ekin_rot() 
    return ekin



def md_step(syst, mol, el, ham, therm, dt, integrator, epot):

    S = MATRIX3x3
    I = MATRIX3x3

    e_t = ekin_tr(syst)
    e_r = ekin_rot(syst)
    therm.update_thermostat_forces(e_t, e_r, 0.0)
    therm.propagate_nhc(0.5*dt, e_t, e_r, 0.0)  # for N-H
    therm.propagate_sPs(0.5*dt)  # for N-P

    s_var = therm.get_s_var()
    dt_half_s = 0.5*dt*s_var
    dt_over_s = (dt/s_var)
    dt_over_s2 = (dt_over_s/s_var)


    S = 0.0
    I.identity()

    S = S + therm.get_ksi_t() * I
    sc1 = (exp_(S,-0.5*dt))
    sc2 = 0.5*dt*(exp1_(S,-0.25*dt))

    #------------------- Angular momentum propagation -----------------------    
    ksi_r = therm.get_ksi_r()
    sc3 = exp(-0.5*dt*ksi_r)
    sc4 = 0.5*dt*exp(-0.25*dt*ksi_r)*sinh_(0.25*dt*ksi_r);


    for n in exrange(syst.Number_of_fragments):
        #-------------------- Linear momentum propagation --------------------
        syst.Fragments[n].Group_RB.scale_linear_(sc1)
        syst.Fragments[n].Group_RB.apply_force(sc2)

        #------------------- Angular momentum propagation -----------------------
        syst.Fragments[n].Group_RB.scale_angular_(sc3)
        syst.Fragments[n].Group_RB.apply_torque(sc4)

    therm.propagate_Ps(-0.5*dt*epot)
    


    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    ps = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Propagate translational DOFs 
        syst.Fragments[n].Group_RB.shift_position(dt * syst.Fragments[n].Group_RB.rb_p * syst.Fragments[n].Group_RB.rb_iM);
 
        # Propagate rotational DOFs
        if integrator=="Jacobi":
            syst.Fragments[n].Group_RB.propagate_exact_rb(dt)
        elif integrator=="DLML":
            ps = syst.Fragments[n].Group_RB.propagate_dlml(dt)

    # Update atomic coordinates
    for n in xrange(syst.Number_of_fragments):
        syst.update_atoms_for_fragment(n)

    #=========== Potential ===================
    syst.zero_forces_and_torques()
    syst.extract_atomic_q(mol.q)

    epot = compute_potential_energy(mol, el, ham, 1)  # 1 - FSSH forces
    compute_forces(mol, el, ham, 1)

    syst.set_atomic_f(mol.f)
    syst.update_fragment_forces_and_torques()


    ekin = 0.0
    for n in xrange(syst.Number_of_fragments):
        # Linear momentum propagation:
        syst.Fragments[n].Group_RB.apply_force(0.5*dt)
        # Angular momentum propagation:
        syst.Fragments[n].Group_RB.apply_torque(0.5*dt)        
        ekin += (syst.Fragments[n].Group_RB.ekin_tr() + syst.Fragments[n].Group_RB.ekin_rot() )

    etot = ekin+epot
    return [ekin, epot, etot]






##############################################################

# Create Universe and populate it
U = Universe()
verbose = 0
Load_PT(U, "elements.dat", verbose)


# Create force field
uff = ForceField({"bond_functional":"Harmonic",
                  "angle_functional":"Fourier",
                  "dihedral_functional":"General0",
                  "oop_functional":"Fourier",
                  "mb_functional":"LJ_Coulomb","R_vdw_on":10.0,"R_vdw_off":15.0 })



Load_UFF(uff)
verb = 0
assign_rings = 1



#======= System ==============
#for i in range(1,13):
for i in [1]:
    print "=================== System ",i,"======================="

    syst = System()
    Load_Molecule(U, syst, os.getcwd()+"/Clusters/23waters.ent", "pdb")
    #Load_Molecule(U, syst, os.getcwd()+"/Clusters/23waters_noq.ent", "pdb")

    
    syst.determine_functional_groups(1)  # 
    syst.show_atoms()

    syst.init_fragments()
#    syst.show_fragments()
#    syst.show_molecules()

    print "Number of atoms in the system = ", syst.Number_of_atoms
    atlst1 = range(1,syst.Number_of_atoms+1)


    # Creating Hamiltonian
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, verb, assign_rings)  

    ham.show_interactions_statistics()

   
    ham.set_system(syst)
    ham.compute()



#    sys.exit(0)

    print "Energy = ", ham.H(0,0), " a.u."
    print "Force 1 = ", ham.dHdq(0,0,0), " a.u."
    print "Force 3 = ", ham.dHdq(0,0,3), " a.u."


#    sys.exit(0)

    #--------------------- Molecular dynamics ----------------------
    el = Electronic(1,0)

    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    print "q= ", mol.q[0]
    print "p= ", mol.p[0]
    print "f= ", mol.f[0]
    print "mass= ", mol.mass[0]
    

    init_md(syst, mol, el, ham)
    syst.init_fragments()
    syst.show_fragments()



    #=================== Propagation ====================

    integrator = "DLML"


    ########################## Cooling #################################

    f = open("_en_cooling.txt","w")
    f.close()
    dt = 1.0 

    for i in xrange(100):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_cooling.xyz",i)

        for j in xrange(10):
            ekin, epot, etot = md_step(syst, mol, el, ham, dt, integrator)

        syst.cool()

        f = open("_en_cooling.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, ekin, epot, etot))
        f.close()


    ########################## Production MD #################################

    # Thermostat
    therm = Thermostat({"Q":100.0, "nu_therm":0.01, "NHC_size":2, "Temperature":300.0, "thermostat_type":"Nose-Hoover" })   
    therm.set_Nf_t(syst.Nf_t); 
    therm.set_Nf_r(syst.Nf_r);
    therm.init_nhc();


    f = open("_en_md.txt","w")
    f.close()
    dt = 20.0 

    syst.init_fragment_velocities(300.0 )

    for i in xrange(100):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)

        for j in xrange(10):
            ekin, epot, etot = md_step_nvt(syst, mol, el, ham, therm, dt, integrator)

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, ekin, epot, etot))
        f.close()




