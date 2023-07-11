#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
###################################################################
# This is a classical all-atomic MD followed by the EHT calculations on each fragment
###################################################################

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

kb = (1.9872065E-3 / 627.5094709)  # in Ha/(mol*K)    

def operator_A(syst, therm, baro, dt_half):

    I = MATRIX3x3(); 
    I.xx = 1.0;  I.yy = 1.0; I.zz = 1.0

    S = baro.ksi_eps_iso * I + (3.0*baro.ksi_eps_iso/baro.get_Nf_t())*I;    
    S = S + therm.get_ksi_t() * I; 

    sc1 = exp_(S,-dt_half);
    sc2 = dt_half*exp1_(S,-dt_half*0.5)

    #print sc1.xx, sc1.xy, sc1.xz, sc1.yx, sc1.yy, sc1.yz, sc1.zx, sc1.zy, sc1.zz
    #print sc2.xx, sc2.xy, sc2.xz, sc2.yx, sc2.yy, sc2.yz, sc2.zx, sc2.zy, sc2.zz

    ksi_r = therm.get_ksi_r();
    sc3 = math.exp(-dt_half*ksi_r);
    sc4 = dt_half*math.exp(-0.5*dt_half*ksi_r)*sinh_(0.5*dt_half*ksi_r);
    #print "sc3 = ", sc3
    #print "sc4 = ", sc4 #sc4.xx, sc4.xy, sc4.xz, sc4.yx, sc4.yy, sc4.yz, sc4.zx, sc4.zy, sc4.zz

    for j in xrange(syst.Number_of_fragments):
        syst.Fragments[j].Group_RB.scale_linear_(sc1)
        syst.Fragments[j].Group_RB.apply_force(sc2)
        syst.Fragments[j].Group_RB.scale_angular_(sc3)
        syst.Fragments[j].Group_RB.apply_torque(sc4)



def operator_B(syst, therm, baro, ST, dt_half):

    baro.update_barostat_forces(syst.ekin_tr(),syst.ekin_rot(),ST.curr_V,ST.curr_P)
    scl = therm.get_ksi_b()
    baro.propagate_velocity(dt_half,scl)


def operator_NHCB(syst, therm, baro, dt_half):
    ekin_baro = baro.ekin_baro()
    therm.update_thermostat_forces(syst.ekin_tr(), syst.ekin_rot(),ekin_baro)
    therm.propagate_nhc(dt_half,syst.ekin_tr(),syst.ekin_rot(),ekin_baro)


def NPT_integrator(ham, el, syst, therm, baro, md, ST, H0):

#    print ">>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<"

    ndof = 3.0 * syst.Number_of_fragments
#    print "ST.curr_V = ", ST.curr_V, " ST.curr_P = ", ST.curr_P
#    print "syst.volume() = ", syst.volume()
#    print baro.get_Nf_b()

    dt_half = 0.5*md.dt

    operator_NHCB(syst, therm, baro, dt_half)
    therm.propagate_sPs(dt_half); 
    operator_B(syst, therm, baro, ST, dt_half)

    s_var = therm.get_s_var(); 
    dt_half_s = dt_half*s_var;
    dt_over_s = (2.0*dt_half/s_var);
    dt_over_s2 = (dt_over_s/s_var);

    operator_A(syst, therm, baro, dt_half)
    therm.propagate_Ps(-dt_half*ST.E_pot)

    sc1 = baro.pos_scale(2.0*dt_half)
    sc2 = 2.0*dt_half*baro.vpos_scale(2.0*dt_half);

    for j in xrange(syst.Number_of_fragments):      
      therm.propagate_Ps( 0.5*dt_over_s2*(syst.Fragments[j].Group_RB.ekin_tr()) )
      s_var = therm.s_var;
      therm.propagate_Ps( 0.5*dt_over_s2*(syst.Fragments[j].Group_RB.ekin_tr()) )
      syst.Fragments[j].Group_RB.scale_position(sc1); 
      syst.Fragments[j].Group_RB.shift_position(sc2 * syst.Fragments[j].Group_RB.rb_iM * syst.Fragments[j].Group_RB.rb_p)

    therm.propagate_Ps(2.0*dt_half*(H0 - ndof * kb * therm.Temperature*(math.log(therm.s_var)+1.0) ) )


#    baro.show_info()
    eps = baro.ksi_eps
#    print eps.xx, eps.xy, eps.xz, eps.yx, eps.yy, eps.yz, eps.zx, eps.zy, eps.zz
    eps = baro.G_eps
#    print eps.xx, eps.xy, eps.xz, eps.yx, eps.yy, eps.yz, eps.zx, eps.zy, eps.zz
#    print baro.eps_iso, baro.ksi_eps_iso, baro.G_eps_iso


    syst.Box = baro.pos_scale(2*dt_half) * syst.Box;

    for j in xrange(syst.Number_of_fragments):
        syst.update_atoms_for_fragment(j)

    syst.zero_forces_and_torques();

    mol = Nuclear(3*syst.Number_of_atoms);
    syst.extract_atomic_q(mol.q); # syst -> mol
    ST.E_pot = compute_forces(mol, el, ham, 1); 
    syst.set_atomic_f(mol.f);    # mol -> syst
    syst.update_fragment_forces_and_torques()

  
#    print "ST.E_pot = ",ST.E_pot
    for j in xrange(syst.Number_of_atoms):    
        x = syst.Atoms[j].Atom_RB.rb_cm.x
        y = syst.Atoms[j].Atom_RB.rb_cm.y
        z = syst.Atoms[j].Atom_RB.rb_cm.z
        fx = syst.Atoms[j].Atom_RB.rb_force.x
        fy = syst.Atoms[j].Atom_RB.rb_force.y
        fz = syst.Atoms[j].Atom_RB.rb_force.z

#        print j, x, y, z, fx, fy, fz



    operator_A(syst, therm, baro, dt_half)

    ST.E_kin = 0.0
    for j in xrange(syst.Number_of_fragments):
        ST.E_kin += syst.Fragments[j].Group_RB.ekin_tr()

    therm.propagate_Ps( -dt_half*ST.E_pot) 


    #=============== Ideal gas virial ==============
    tmp = MATRIX3x3()
    Pid = MATRIX3x3()
    p_tot = VECTOR()
    m_tot = 0.0
    for j in xrange(syst.Number_of_atoms):
        pj = syst.Atoms[j].Atom_RB.rb_p
        tmp.tensor_product(pj, pj);
        Pid = Pid + syst.Atoms[j].Atom_RB.rb_iM * tmp;
        p_tot = p_tot + pj; 
        m_tot = m_tot + syst.Atoms[j].Atom_RB.rb_mass;

    tmp.tensor_product(p_tot,p_tot) 
    Pid = Pid - (1.0/m_tot) * tmp;

    print "Pid = "
    print Pid.xx, Pid.xy, Pid.xz
    print Pid.yx, Pid.yy, Pid.yz
    print Pid.zx, Pid.zy, Pid.zz
   
    Pvir = ham.get_stress("at") 
    print "Virial = "
    print Pvir.xx, Pvir.xy, Pvir.xz
    print Pvir.yx, Pvir.yy, Pvir.yz
    print Pvir.zx, Pvir.zy, Pvir.zz

    X  =  (Pid + Pvir) /syst.volume() # total pressure tensor!

    ST.curr_P_tens = 0.5* (X + X.T() )

    ST.curr_P = (ST.curr_P_tens.xx + ST.curr_P_tens.yy + ST.curr_P_tens.zz)/3.0    
    ST.curr_V = syst.volume();
    print "ST.curr_P_tens = "
    print ST.curr_P_tens.xx, ST.curr_P_tens.xy, ST.curr_P_tens.xz
    print ST.curr_P_tens.yx, ST.curr_P_tens.yy, ST.curr_P_tens.yz
    print ST.curr_P_tens.zx, ST.curr_P_tens.zy, ST.curr_P_tens.zz
    print "ST.curr_P = ", ST.curr_P
    print "ST.curr_V = ", ST.curr_V

    operator_B(syst, therm, baro, ST, dt_half)
    therm.propagate_sPs(dt_half)
    operator_NHCB(syst, therm, baro, dt_half)


    ST.E_kin = ST.E_kin/(therm.s_var*therm.s_var)
    ST.E_kin_tr = syst.ekin_tr()
    ST.E_kin_rot = syst.ekin_rot()
    ST.E_tot = ST.E_kin + ST.E_pot
    ST.H_NP = ST.E_tot + baro.ekin_baro() + ST.curr_V * baro.Pressure
    ST.curr_T = 2.0*ST.E_kin/(ndof* kb)



def main():

    rnd = Random()

    #--------------------- Initialization ----------------------

    # Create Universe and populate it
    U = Universe(); LoadPT.Load_PT(U, os.getcwd()+"/elements.txt")

    # Create force field
#    uff = ForceField({"mb_functional":"LJ_Coulomb","R_vdw_on":40.0,"R_vdw_off":55.0 })  # this can not be used for PBC
    uff = ForceField({"R_vdw_on":10.0,"R_vdw_off":12.0, "mb_functional":"vdw_LJ1","mb_excl_functional":"vdw_LJ1"})

    LoadUFF.Load_UFF(uff,"uff.dat")


    # Create molecular system and initialize the properties
    syst = System()

    LoadMolecule.Load_Molecule(U, syst, "au.pdb", "true_pdb")



    syst.determine_functional_groups(0)  # do not assign rings
    syst.init_fragments()
    print "Number of atoms in the system = ", syst.Number_of_atoms
    print "Number of bonds in the system = ", syst.Number_of_bonds
    print "Number of angles in the system = ", syst.Number_of_angles
    print "Number of dihedrals in the system = ", syst.Number_of_dihedrals
    print "Number of impropers in the system = ", syst.Number_of_impropers
    atlst1 = range(1,syst.Number_of_atoms+1)

    T1 =  VECTOR(32.6970772436, 0.0, 0.0)
    T2 =  VECTOR(16.3485386218, 28.3164995224, 0.0)
    T3 =  VECTOR(0.0, 0.0, 26.6970517757)
    syst.init_box(T1, T2, T3)

    # Creating Hamiltonian and initialize it
    ham = Hamiltonian_Atomistic(1, 3*syst.Number_of_atoms)
    ham.set_Hamiltonian_type("MM")
    ham.set_interactions_for_atoms(syst, atlst1, atlst1, uff, 1, 0)  # 0 - verb, 0 - assign_rings
    ham.show_interactions_statistics()


    # Bind Hamiltonian and the system   
    ham.set_system(syst);   ham.compute();   print "Energy = ", ham.H(0,0), " a.u."

    # Electronic DOFs
    el = Electronic(1,0)


    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)

    # Initialize MD variables
    nve_md.nve_md_init(syst, mol, el, ham)



    #=================== Propagation ====================
 
    ########################## Cooling #################################

    md = MD({"max_step":1,"ensemble":"NPT","integrator":"DLML","terec_exp_size":10,"dt":20.0,"n_medium":1,"n_fast":1,"n_outer":1})
    md.show_info()
    
    # Thermostat
    therm = Thermostat({"Temperature":278.0,"Q":100.0,"thermostat_type":"Nose-Hoover","nu_therm":0.001,"NHC_size":1})
    therm.show_info()

    # Barostat
    baro = Barostat({"W":10000.0,"Pressure":1.0,"nu_baro":0.001})
    baro.show_info()



    ST = State() 
    ST.set_system(syst)
    ST.set_thermostat(therm)
    ST.set_barostat(baro)
    ST.set_md(md)

    ST.init_md(mol, el, ham, rnd)    



    H0 = ST.E_tot; 
    H0 = H0 + therm.energy()

    f = open("_en_cooling.txt","w")
    f.close()

    md.max_step = 50
    for i in xrange(100):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_cooling.xyz",i)

        #ST.run_md(mol, el, ham)
#        print i
#        ST.run_md(el, ham)
        NPT_integrator(ham, el, syst, therm, baro, md, ST, H0)
        ekin = ST.E_kin;    epot = ST.E_pot
        ST.cool()

        f = open("_en_cooling.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f  curr_P= %10.8f  curr_V= %8.5f\n" 
                % (i, ST.E_kin, ST.E_pot, ST.E_tot, ST.H_NP, ST.curr_T, ST.curr_P, ST.curr_V ))
        f.close()

#    sys.exit(0)

    ########################## Production MD #################################

    syst.init_atom_velocities(300.0, rnd)  # must be this !!!
#    syst.init_fragment_velocities(300.0, rnd)

    f = open("_en_md.txt","w")
    f.close()
    md.dt = 40.0 
    md.max_step = 10

#    H0 = ST.E_tot; 
#    H0 = H0 + therm.energy()

    print ">>>>>>>>>>>>> <<<<<<<<<<<<<"

    for i in xrange(500):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)

        NPT_integrator(ham, el, syst, therm, baro, md, ST, H0)
#        ST.run_md(mol, el, ham)

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f  curr_P= %10.8f  curr_V= %8.5f\n" 
                % (i, ST.E_kin, ST.E_pot, ST.E_tot, ST.H_NP, ST.curr_T, ST.curr_P, ST.curr_V ))
        f.close()


main()

