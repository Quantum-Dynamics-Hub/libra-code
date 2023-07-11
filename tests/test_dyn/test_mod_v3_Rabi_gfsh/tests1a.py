import os
import sys
import math
import cmath

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



rnd = Random()



def boltz():
# This function generates x and p taken from the Boltzmann distribution
    #rnd = Random() - must call it only once outside!!!

    kT = 9.5e-4         # 300.0 K or so
    p_ = math.sqrt(kT)  # momentum that corresponds to temperatures
    
    Eold = 0.5*p_*p_    # energy, m = 1 - electron
    Enew = 0.0


    for i in xrange(1000):
        p = p_ + 1.1*math.sqrt(kT)*rnd.uniform(-1.0, 1.0) #rnd.normal() # proposed state p

        Enew = 0.5*p*p
        dE = Enew - Eold

        ksi = rnd.uniform(0.0,1.0)
        prob = 1.0
        argg = dE/kT
        if argg >40:
            prob = 0.0
        elif argg < -40:
            prob = 1.0
        else: 
            prob = math.exp(-dE/kT)  #min(1.0, math.exp(-dE/kT))

        if(ksi<prob):  # accept new state with Metropolis scheme
            Eold = Enew
            p_ = p

    return p_




def q_ave(ntraj, mol, indx):

    x = 0.0
    for j in xrange(ntraj): # for all trajectories
        x = x + mol[j].q[indx]
    x = x / float(ntraj)

    sigma = 0
    for j in xrange(ntraj): # for all trajectories
        sigma = sigma + (mol[j].q[indx] - x)**2

    sigma = math.sqrt(sigma / float(ntraj))

    
    return x, sigma

def p_ave(ntraj, mol, indx):

    p = 0.0
    for j in xrange(ntraj): # for all trajectories
        p = p + mol[j].p[indx]
    p = p / float(ntraj)

    sigma = 0
    for j in xrange(ntraj): # for all trajectories
        sigma = sigma + (mol[j].p[indx] - p)**2

    sigma = math.sqrt(sigma / float(ntraj))
    
    return p, sigma




def sh_pop1(ntraj, nelec, el):

    pops = [0.0, 0.0]
 
    for j in xrange(ntraj): # for all trajectories
        pops[ el[j].istate ] += 1

    for i in xrange(nelec):
        pops[i] = pops[i] / float(ntraj)

    return pops


def se_pop1(ntraj, nelec, el):

    pops = [0.0, 0.0, 0.0, 0.0]
 
    for j in xrange(ntraj): # for all trajectories
        pops[0] = pops[0] + el[j].rho(0,0).real
        pops[1] = pops[1] + el[j].rho(1,1).real
        pops[2] = pops[2] + el[j].rho(0,1).real
        pops[3] = pops[3] + el[j].rho(0,1).imag

    for i in xrange(4):
        pops[i] = pops[i] / float(ntraj)


    return pops




def apply_thermostat(MOL, THERM, therms, S, ntraj, t_method, dt):

    for i in xrange(ntraj):
        if t_method=="nhc":
            MOL[i].p[0] = MOL[i].p[0] * THERM.vel_scale(dt)
        elif t_method=="ens_nhc":
            MOL[i].p[0] = MOL[i].p[0] * therms[i].vel_scale(dt)
        elif t_method=="lang":
            MOL[i].p[0] = MOL[i].p[0] * math.exp(-gamma*0.5*dt)
            MOL[i].p[0] = MOL[i].p[0] + sigma * S[i] * dt
            MOL[i].p[0] = MOL[i].p[0] * math.exp(-gamma*0.5*dt)




print "\nTest 2: Set parameters"
#######################################
# 1)
#method = "fssh"
method = "gfsh"
#method = "esh"
#method = "mssh"
# 2)
#t_method = "lang" 
#t_method = "nhc"
#t_method = "ens_nhc"
t_method = "none"
# 3)
do_reverse = 0
# 4)
nu_therm = 0.01
# 5)
#rang = range(0,17)
#rang = range(10,17)
rang = [0]

#6) 
sh_opt = 1  # 0 - MF, 1 - SH
ntraj = 2500 # 25000

#######################################

E0 = -0.001
E1 = 0.001
V01= 0.001 # * math.sqrt(10.0)
D = 0.0019
Lx = 1.0
Ly = 1.0

nsnap = 200
nstep = 250
dt = 0.1
isurface = 1
rep = 0 # diabatic
ham_indx = 200 # Rabi modified



use_boltz_factor = 0
T = 300.0
kb = 3.166811429e-6
do_rescaling = 1 
nel = 2
nnucl = 2


mass = 2000.0
omega = 3.5e-4
Er = 2.39e-2
V_ = 5.0e-5
kT = 9.5e-4
dE = (3.0e-2 - 1.5e-2)/16.0
gamma = 3.0e-4
#gamma = 3.0e-5
sigma = math.sqrt(2.0*mass*gamma*kT/dt)




A = 0.5/math.sqrt(1.0 + 0.25*((E1-E0)/V01)**2)
OMEGA = math.sqrt(V01**2 + 0.25*(E1-E0)**2)


#f = math.exp(-math.sqrt((E1-E0)**2 + 4.0*V01**2)/(kb*T))
f = math.exp(-(E1-E0)/(kb*T))
print "Boltzmann ratio: ", f, "p0(eq)= ", 1.0/(1.0 + f), "p1(eq)= ", f/(1.0+f)

for ie in rang: 

    initstate = intList()
    finstate = intList()
    glist = MATRIXList()
    mollist = NuclearList()
    hamlist = HamiltonianList()
    ksilist = doubleList()

    for i in xrange(ntraj):
        initstate.append(0)
        finstate.append(0)
        g = MATRIX(2,2)
        glist.append(g)
        ksilist.append(0.0)


    print "ie = ", ie
    print "Step 1: Initialize ensemble and Hamiltonians"
 
    #----------------------- Initialization ---------------------------
    # Create Electronic DOFs
    EL = []
    for i in xrange(ntraj):
        el = Electronic(nel,0)
        EL.append(el)

    # Create Nuclear DOFs
    MOL = []
    for i in xrange(ntraj):
        mol = Nuclear(nnucl)
        mol.mass[0] = mass
        mol.mass[1] = mass
        mol.q[0], mol.p[0], mol.f[0] = 0.0, 2.0+boltz(), 0.0
        mol.q[1], mol.p[1], mol.f[1] = 0.0, 2.0+boltz(), 0.0
        MOL.append(mol)
        mollist.append(mol)


    # Create Hamiltonians
    HAM = []
    for i in xrange(ntraj):
        ham = Hamiltonian_Model(ham_indx)  
        ham.set_rep(rep)
        ham.set_params([E0, E1, V01, D, Lx, Ly ])
        ham.set_v( [ MOL[i].p[0]/MOL[i].mass[0], MOL[i].p[1]/MOL[i].mass[1] ])
        ham.set_q( [ MOL[i].q[0], MOL[i].q[1] ] )
        HAM.append(ham) 
        hamlist.append(ham)


    # Create thermostats 
    # The entangled thermostat!
    THERM = Thermostat({"nu_therm":nu_therm, "NHC_size":15, "Temperature":300.0, "thermostat_type":"Nose-Hoover" })   
    THERM.set_Nf_t(1)
    THERM.set_Nf_r(0);
    THERM.init_nhc();


    # An ensemble of disentangled thermostats
    therms = []  
    for i in xrange(ntraj):
        therm = Thermostat({"nu_therm":nu_therm, "NHC_size":2, "Temperature":300.0, "thermostat_type":"Nose-Hoover" }) 
        therm.set_Nf_t(1);
        therm.set_Nf_r(0);
        therm.init_nhc();
        therms.append(therm)


    # Initialize electronic states - this depends on the state of Hamiltonian
    for i in xrange(ntraj):   
        EL[i].istate = isurface

        if isurface==0:
            f0, f1 = 1.0, 0.0
        else:
            f0, f1 = 0.0, 1.0
        ksi = rnd.uniform(-2.0*math.pi, 2.0*math.pi);
        EL[i].q[0] = f0 * math.cos(ksi)
        EL[i].p[0] = f0 * math.sin(ksi)
        EL[i].q[1] = f1 * math.cos(ksi)
        EL[i].p[1] = f1 * math.sin(ksi)



    # Update energies and forces for this setup
    e_tot = 0.0
    e_pot = 0.0
    e_kin = 0.0
    EPOT = []
    EKIN = []
    EBAT = []
    for i in xrange(ntraj):
#        epot = compute_potential_energy(MOL[i], EL[i], HAM[i], 1)  # 0 - MF forces, 1 - FSSH (state-resolved) forces
        epot = compute_forces(MOL[i], EL[i], HAM[i], sh_opt) # 0 - MF, 1 - FSSH
        ekin = compute_kinetic_energy(MOL[i])
        EPOT.append(epot)
        EKIN.append(ekin)
        e_pot = e_pot + epot
        e_kin = e_kin + ekin

    e_pot = e_pot / float(ntraj)
    e_kin = e_kin / float(ntraj)
    e_tot = e_pot + e_kin
    
    e_bath = 0.0
    if t_method=="nhc":
        e_bath = THERM.energy()
    elif t_method=="ens_nhc":
        for i in xrange(ntraj):
            e_bath = e_bath + therms[i].energy()
        e_bath = e_bath / float(ntraj)


    e_ext = e_tot + e_bath



    t = 0.0
    pop0_ex = (2.0*A*math.sin(OMEGA*t))**2
    pops = sh_pop1(ntraj, nel, EL)
    pops_se = se_pop1(ntraj, nel, EL)
#    x, sigmax_ = x_ave(ntraj, MOL)
#    p, sigmap_ = p_ave(ntraj, MOL)

    x, sigmax_ = q_ave(ntraj, MOL, 0)
    px, sigmapx_ = p_ave(ntraj, MOL, 0)
    y, sigmay_ = q_ave(ntraj, MOL, 1)
    py, sigmapy_ = p_ave(ntraj, MOL, 1)



    snap = 0
    f = open("relax"+str(ie)+".txt","w")
#    f.write("%10.8f   %10.8f   %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n" % (snap*nstep*dt, pops[0], pops[1], x, sigmax_, e_kin, e_pot, e_tot, e_ext, pops_se[0], pops_se[1], p, sigmap_, pops_se[2], pops_se[3], pop0_ex))    
    f.write("%10.8f   %10.8f   %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n" % (t, pops[0], pops[1], x, y, sigmax_, sigmay_, e_kin, e_pot, e_tot, e_ext, pops_se[0], pops_se[1], px, py, sigmapx_, sigmapy_, pops_se[2], pops_se[3], pop0_ex))    
    f.close()


    #---------------------- Propagation -------------------------------    
    S = []
    for i in xrange(ntraj):
        S.append(0.0)

    print "Step 2: Propagation"
    for snap in xrange(nsnap):
        for step in xrange(nstep):

            for i in xrange(ntraj):
                EL[i].propagate_electronic(0.5*dt,HAM[i])

        
            #******************** Nuclear dynamics *****************************
            #----- Thermostat -------
            for i in xrange(ntraj):
                S[i] = rnd.normal()

            if t_method=="nhc" or t_method=="ens_nhc":
                apply_thermostat(MOL, THERM, therms, S, ntraj, t_method, 0.5*dt)



            e_pot = 0.0
            e_kin = 0.0
            for i in xrange(ntraj):
                #============= Propagate nuclear DOFs =============
                MOL[i].propagate_p(0.5*dt)
                MOL[i].propagate_q(dt)

                #========= Update forces and energies =============
                #EPOT[i] = compute_potential_energy(MOL[i], EL[i], HAM[i], 1)  # 0 - MF forces, 1 - FSSH (state-resolved) forces
                EPOT[i] = compute_forces(MOL[i], EL[i], HAM[i], sh_opt) # 0 - MF, 1 - FSSH
                e_pot = e_pot + EPOT[i]

                EKIN[i] = compute_kinetic_energy(MOL[i])
                e_kin = e_kin + EKIN[i]

            e_pot = e_pot / float(ntraj)
            e_kin = e_kin / float(ntraj)


            #========== Propagate thermostat ==================
            if t_method=="nhc":
                THERM.propagate_nhc(dt, e_kin, 0.0, 0.0)
            elif t_method=="ens_nhc":
                for i in xrange(ntraj):                    
                    therms[i].propagate_nhc(dt, EKIN[i], 0.0, 0.0)

            
            for i in xrange(ntraj):
                #============= Propagate nuclear DOFs =============
                MOL[i].propagate_p(0.5*dt)


            if t_method=="nhc" or t_method=="ens_nhc":
                apply_thermostat(MOL, THERM, therms, S, ntraj, t_method, 0.5*dt)


            #******************** End of nuclear dynamics *****************************

            e_kin = 0.0
            for i in xrange(ntraj):
                HAM[i].set_v([MOL[i].p[0]/MOL[i].mass[0] ])
                EL[i].propagate_electronic(0.5*dt,HAM[i])  
                EKIN[i] = compute_kinetic_energy(MOL[i])
                e_kin = e_kin + EKIN[i]
            e_kin = e_kin / float(ntraj)
            e_tot = e_pot + e_kin




        
            #***************** Now, incorporate surface hop  ****************************
            if sh_opt==1:

                g = MATRIX(2,2)
                
                e_pot = 0.0        
                for i in xrange(ntraj):
                    if method=="mssh":
                        compute_hopping_probabilities_mssh(MOL[i], EL[i], HAM[i], g, dt, use_boltz_factor, T)
                    elif method=="fssh":
                        compute_hopping_probabilities_fssh(MOL[i], EL[i], HAM[i], g, dt, use_boltz_factor, T)
                    elif method=="gfsh":
                        compute_hopping_probabilities_gfsh(MOL[i], EL[i], HAM[i], g, dt, use_boltz_factor, T)
                
                    ksi = rnd.uniform(0.0, 1.0)
                
                    ksilist[i] = ksi
                    initstate[i] = EL[i].istate 
                    mollist[i] = MOL[i]
                    hamlist[i] = HAM[i]
                    glist[i] = g
                
                    if method=="mssh" or method=="fssh" or method=="gfsh":
                        old = EL[i].istate
                        #print "Momentum before hop= ", MOL[i].p[0]
                        EL[i].istate = hop(EL[i].istate, MOL[i], HAM[i], ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary
                        #print "Momentum after hop= ", MOL[i].p[0]
                
                
                if method=="esh":
                    finstate = hop(ntraj, initstate, mollist, hamlist, ksilist, glist, do_rescaling, rep, do_reverse)
                    for i in xrange(ntraj):
                        EL[i].istate = finstate[i]
                        MOL[i] = mollist[i]
                
                
                for i in xrange(ntraj):
                    # Hop has happened - recompute energy and forces    
                    #if EL[i].istate != old:
                    #EPOT[i] = compute_potential_energy(MOL[i], EL[i], HAM[i], 1)  # 0 - MF forces, 1 - FSSH (state-resolved) forces
                    EPOT[i] = compute_forces(MOL[i], EL[i], HAM[i], sh_opt) # 0 - MF, 1 - FSSH
                    e_pot = e_pot + EPOT[i]  # Either use old values, or the updated ones
                e_pot = e_pot / float(ntraj)
                e_tot = e_pot + e_kin



        
        #---------------------- Compute and print info for given snap ------------------------------
        t = (snap+1)*nstep*dt
        pop0_ex = (2.0*A*math.sin(OMEGA*t))**2

        pops = sh_pop1(ntraj, nel, EL)
        pops_se = se_pop1(ntraj, nel, EL)
        x, sigmax_ = q_ave(ntraj, MOL, 0)
        px, sigmapx_ = p_ave(ntraj, MOL, 0)
        y, sigmay_ = q_ave(ntraj, MOL, 1)
        py, sigmapy_ = p_ave(ntraj, MOL, 1)



        e_bath = 0.0
        if t_method=="nhc":
            e_bath = THERM.energy()
        elif t_method=="ens_nhc":
            for i in xrange(ntraj):
                e_bath = e_bath + therms[i].energy()
            e_bath = e_bath / float(ntraj)

        e_ext = e_tot + e_bath


        f = open("relax"+str(ie)+".txt","a")
        f.write("%10.8f   %10.8f   %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n" % (t, pops[0], pops[1], x, y, sigmax_, sigmay_, e_kin, e_pot, e_tot, e_ext, pops_se[0], pops_se[1], px, py, sigmapx_, sigmapy_, pops_se[2], pops_se[3], pop0_ex))    

#        f.write("%10.8f   %10.8f   %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  
#        %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n" %
#        (t, pops[0], pops[1], x, y, sigmax_, sigmay_, e_kin, e_pot, e_tot, e_ext, pops_se[0], pops_se[1],
#        px, py, sigmapx_, sigmapy_, pops_se[2], pops_se[3], pop0_ex))    

        f.close()

    

print "time propagation is done"      


