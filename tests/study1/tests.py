import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygdyn import *
from cyghamiltonian import *



print "\nTest 2: Set parameters"
method = "esh"
ntraj = 50
use_boltz_factor = 0
T = 300.0
kb = 3.166811429e-6
rnd = Random()
do_rescaling = 1 
do_reverse = 1
rep = 1 # adiabatic
ham_indx = 2
nel = 2
nnucl = 1
dt = 1.0



P_list = []
for i in xrange(40):
    P_list.append(1.0 + 1.0*i)



f = open("p_scatt_re_"+method+".txt","w")
f.close()
f = open("p_scatt_tr_"+method+".txt","w")
f.close()


for P in P_list:

    dt = 1.0 #0.01*(2000.0/P)

    T =  0.5*P*P/(2000.0*kb)

    print "P = ", P
    print "Step 1: Initialize ensemble and Hamiltonians"
    ens = Ensemble(ntraj, nel, nnucl) 
    ens.ham_set_ham("model", ham_indx)
    ens.ham_set_rep(rep)

    # Initialization
    for i in xrange(ntraj):
        ens.mol[i].mass[0] = 2000.0
        ens.mol[i].q[0] = -14.0
        ens.mol[i].p[0] = P + 0.01*(rnd.normal() - 0.5)*P
        ens.ham_set_v(i, [ens.mol[i].p[0]/ens.mol[i].mass[0]] )

    epot = compute_potential_energy(ens, 1)  # 0 - MF forces, 1 - FSSH
    compute_forces(ens, 1)

    #---------------------- Propagation -------------------------------
    print "Step 2: Propagation"
    for t in xrange(5000):
        ens.el_propagate_electronic(0.5*dt)

        # Nuclear dynamics
        ens.mol_propagate_p(0.5*dt)
        ens.mol_propagate_q(dt)
        ens.ham_set_v()

        epot = compute_potential_energy(ens, 1)  # 0 - MF forces, 1 - FSSH
        compute_forces(ens, 1)

        ens.mol_propagate_p(0.5*dt)

        # el-dyn
        ens.el_propagate_electronic(0.5*dt)    
        ens.ham_set_v()

        ekin = compute_kinetic_energy(ens)


        # Now, incorporate surface hop
        g = MATRIX(2,2)

        if method=="esh":
            compute_hopping_probabilities_esh(ens, g, dt, use_boltz_factor, T)

            for i in xrange(ntraj):
                ksi = rnd.uniform(0.0, 1.0)
                ens.el[i].istate = hop(ens.el[i].istate, ens, i , ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary

        else:  
            for i in xrange(ntraj):
                if method=="mssh":
                    compute_hopping_probabilities_mssh(ens, i, g, dt, use_boltz_factor, T)
                elif method=="fssh":
                    compute_hopping_probabilities_fssh(ens, i, g, dt, use_boltz_factor, T)
                elif method=="gfsh":
                    compute_hopping_probabilities_gfsh(ens, i, g, dt, use_boltz_factor, T) 

                ksi = rnd.uniform(0.0, 1.0)
                ens.el[i].istate = hop(ens.el[i].istate, ens, i , ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary

    #---------------------- Propagation is done ------------------------------

    print "Step 3: Analysis and print"
    se_pops_tr = ens.se_pop(0.0, 100000.0)
    sh_pops_tr = ens.sh_pop(0.0, 100000.0)
    se_pops_re = ens.se_pop(-100000.0, -0.0)
    sh_pops_re = ens.sh_pop(-100000.0, -0.0)


    nrm_se = (1.0/(se_pops_re[0] + se_pops_re[1] + se_pops_tr[0] + se_pops_tr[1]))
    nrm_sh = (1.0/(sh_pops_re[0] + sh_pops_re[1] + sh_pops_tr[0] + sh_pops_tr[1]))

    f = open("p_scatt_re_"+method+".txt","a")
    f.write("P= %8.5f  se_pop_re(0)= %8.5f  se_pop_re(1)= %8.5f  sh_pop_re(0)= %8.5f sh_pop_re(1)= %8.5f\n" % (P, nrm_se*se_pops_re[0], nrm_se*se_pops_re[1], nrm_sh*sh_pops_re[0], nrm_sh*sh_pops_re[1] ))
    f.close()

    f = open("p_scatt_tr_"+method+".txt","a")
    f.write("P= %8.5f  se_pop_tr(0)= %8.5f  se_pop_tr(1)= %8.5f  sh_pop_tr(0)= %8.5f sh_pop_tr(1)= %8.5f\n" % (P, nrm_se*se_pops_tr[0], nrm_se*se_pops_tr[1], nrm_sh*sh_pops_tr[0], nrm_sh*sh_pops_tr[1] ))
    f.close()


print "time propagation is done"      
print "file is closed"


