#*********************************************************************************
#* Copyright (C) 2016-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *
from common import compute_statistics, compute_statistics2, lattice, lattice_p, compute_rmsd, compute_cv
kb = (1.9872065E-3 / 627.5094709)  # in Ha/(mol*K)

def main(params):

    Q, P, ETOT = [], [], []

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

    md = MD({"max_step":10,"ensemble":"NVT","integrator":"DLML","terec_exp_size":10,"dt":20.0,"n_medium":1,"n_fast":1,"n_outer":1})
    md.show_info()
    
    # Thermostat
    therm = Thermostat({"Temperature":params["T"],"Q":100.0,"thermostat_type":"Nose-Hoover","nu_therm":0.001,"NHC_size":5})
    therm.show_info()


    ST = State() 
    ST.set_system(syst)
    ST.set_thermostat(therm)
    ST.set_md(md)

    ST.init_md(mol, el, ham, rnd)    


    f = open("_en_cooling.txt","w")
    f.close()

    md.dt = params["dt_cooling"]
    md.max_step = params["nsteps_cooling"]
    for i in xrange(params["ncycles_cooling"]):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_cooling.xyz",i)

        ST.run_md(mol, el, ham)
        ekin = ST.E_kin;    epot = ST.E_pot
        ST.cool()

        f = open("_en_cooling.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ekin, epot, ST.E_tot, ST.H_NP, ST.curr_T ))
        f.close()


    ########################## Production MD #################################

    syst.init_atom_velocities(params["T"], rnd) 

    f = open("_en_md.txt","w")
    f.close()

    nat = syst.Number_of_atoms
    md.dt = params["dt_md"]
    md.max_step = params["nsteps_md"]
    for i in xrange(params["ncycles_md"]):
        syst.set_atomic_q(mol.q)
        syst.print_xyz("_mol_md.xyz",i)

        ST.run_md(mol, el, ham)

        #==== Sample point =====
        q = MATRIX(3*nat, 1)
        p = MATRIX(3*nat, 1)
        for at in xrange(nat):
            q.set(3*at,   0, mol.q[3*at]);    p.set(3*at,   0, mol.p[3*at]/mol.mass[at])
            q.set(3*at+1, 0, mol.q[3*at+1]);  p.set(3*at+1, 0, mol.p[3*at+1]/mol.mass[at])
            q.set(3*at+2, 0, mol.q[3*at+2]);  p.set(3*at+2, 0, mol.p[3*at+2]/mol.mass[at])
        Q.append( q )
        P.append( p )
        ETOT.append([ST.E_tot])
                 

        f = open("_en_md.txt","a")
        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f  H_NP= %8.5f  curr_T= %8.5f\n" % (i, ST.E_kin, ST.E_pot, ST.E_tot, ST.H_NP, ST.curr_T ))
        f.close()

    return Q, P, ETOT


#==========================================================
params = {}
params["T"] = 278.0

params["dt_cooling"] = 20.0
params["ncycles_cooling"] = 50
params["nsteps_cooling"] = 50

params["dt_md"] = 40.0
params["ncycles_md"] = 11100
params["nsteps_md"] = 10

params["start_step"] = 1100


Q, P, Etot = main(params)

compute_rmsd(Q, params["nsteps_md"]*params["dt_md"], "_rmsd.txt")


cv = compute_cv(Etot, params["start_step"], params["nsteps_md"]*params["dt_md"], "_cv.txt")
print "C_v = ", cv / (kb * params["T"] * params["T"])


vel = []
for p in P:
    vel.append(p.col(0))
acf_matrix.recipe1(vel, params["nsteps_md"]*params["dt_md"]*0.02419, 5000.0, 1.0, "_acf_vel.txt", "_spectrum_vel.txt", 0)




