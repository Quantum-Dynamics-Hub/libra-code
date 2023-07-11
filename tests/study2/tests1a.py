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


def boltz_old():
# This function generates x and p taken from the Boltzmann distribution
    rnd = Random()

    Er = 2.39e-2
    omega = 3.5e-4
    kT = 9.5e-4

 
    mo2 = 0.5*omega*omega # mass = 1
    M = math.sqrt(mo2*Er)

    x_ = -M/mo2           # minimum
    p_ = math.sqrt(kT)         # momentum that corresponds to temperatures

    Eold = mo2*x_*x_ + M*x_ + 0.5*p_*p_    # energy
    Enew = 0.0

    for i in xrange(1000):
        x = 3.0*(M/mo2)*rnd.normal()  # proposed state x
        p = 3.0*math.sqrt(kT)*rnd.normal() # proposed state p

        Enew = mo2*x*x + M*x + 0.5*p*p
        dE = Enew - Eold

        ksi = rnd.uniform(0.0,1.0)
        prob = min(1.0, math.exp(-dE/kT))

        if(ksi<prob):  # accept new state with Metropolis scheme
            Eold = Enew
            x_ = x
            p_ = p

    return [x_, p_]



def boltz():
# This function generates x and p taken from the Boltzmann distribution
    #rnd = Random() - must call it only once outside!!!

    Er = 2.39e-2
    omega = 3.5e-4
    kT = 9.5e-4
    mo2 = 0.5*omega*omega # mass = 1
    M = math.sqrt(mo2*Er)


    X_ = -0.5*M/mo2             # minimum
    p_ = math.sqrt(kT)          # momentum that corresponds to temperatures
    x_ = X_ + 50.0*rnd.normal() # + 0.5*(M/mo2)*rnd.normal()  #(M/mo2)*

    
    Eold = mo2*x_*x_ + M*x_ + 0.5*p_*p_    # energy
    Enew = 0.0


    for i in xrange(1000):
        #print i, x_
        x = x_ + 10.0*rnd.uniform(-1.0, 1.0) #50.0*rnd.normal()  #0.1*(M/mo2)*rnd.normal()  # proposed state x
        p = p_ + 1.1*math.sqrt(kT)*rnd.uniform(-1.0, 1.0) #rnd.normal() # proposed state p

        Enew = mo2*x*x + M*x + 0.5*p*p
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
            x_ = x
            p_ = p

    return [x_, p_]




def x_ave(ens):

    x = 0.0
    for j in xrange(ens.ntraj): # for all trajectories
        x = x + ens.mol[j].q[0]

    x = x / float(ens.ntraj)

    sigma = 0
    for j in xrange(ens.ntraj): # for all trajectories
        sigma = sigma + (ens.mol[j].q[0] - x)**2

    sigma = math.sqrt(sigma / float(ens.ntraj))

    
    return x, sigma



def sh_pop1_ext(ens, rep, xmin, xmax):
# ens - ensemble

    pops = [0.0, 0.0]
 
    for j in xrange(ens.ntraj): # for all trajectories

        ens.ham_set_q(j, ens.mol[j].q )
        #ens.ham_compute_adiabatic(i)
        ens.ham_compute(j)

        # Adiabatic energies
        ens.ham_set_rep(j,1)
        E0 = ens.ham_H(j,0,0).real
        E1 = ens.ham_H(j,1,1).real

        # Extract diabatic energies:
        ens.ham_set_rep(j,0)
        H0 = ens.ham_H(j,0,0).real 
        H1 = ens.ham_H(j,1,1).real
        V  = ens.ham_H(j,0,1).real
        ens.ham_set_rep(j, rep) # set representation to the original one


        for i in xrange(ens.nelec):  # for all electronic states

            if(ens.mol[j].q[0] > xmin and ens.mol[j].q[0] < xmax):  # only those trajectories that are in given range

                if(ens.el[j].istate == 0): # we are in 0-th adiabatic state
        
                    if(i==0):             # Probability to be on 0-th (left) diabat - f0           
                        pops[0] = pops[0] + V*V/((H0 - E0)*(H0 - E0) + V*V);  # f0^2

                    elif(i==1):           # Probability to be on 1-th (right) diabat - g0           
                        pops[1] = pops[1] + (H0 - E0)*(H0 - E0)/((H0 - E0)*(H0 - E0) + V*V);  # g0^2
                # state == 0

                elif(ens.el[j].istate == 1):  # on 1-st adiabatic state

                    if(i==0):             # Probability to be on 0-th (left) diabat - f1         
                        pops[0] = pops[0] + V*V/((H0 - E1)*(H0 - E1) + V*V);  # f1

                    elif(i==1):           # Probability to be on 1-th (right) diabat - g1
                        pops[1] = pops[1] + (H0 - E1)*(H0 - E1)/((H0 - E1)*(H0 - E1) + V*V);  # g1

                # state == 1
            # if  xmin< x <xmax
        # for i
    # for j

    for i in xrange(ens.nelec):  # all electronic states
        pops[i] = pops[i] / float(ens.ntraj)


    return pops




print "\nTest 2: Set parameters"
method = "fssh"
#method = "esh"

t_method = "lang" 
#t_method = "nhc"
#t_method = "ens_nhc"

do_reverse = 0



ntraj = 25
use_boltz_factor = 0
T = 300.0
kb = 3.166811429e-6
do_rescaling = 1 
rep = 1 # adiabatic
ham_indx = 3
nel = 2
nnucl = 1
nsnap = 1000
nstep = 200
dt = 2.5
isurface = 0


mass = 1.0
omega = 3.5e-4
Er = 2.39e-2
V_ = 5.0e-5
kT = 9.5e-4
dE = (3.0e-2 - 1.5e-2)/16.0
#eps_0 = 0.0

gamma = 3.0e-4
sigma = math.sqrt(2.0*mass*gamma*kT/dt)



e_tot = 0.0

for ie in range(0,1):
#for ie in [8]:

    print "ie = ", ie
    print "Step 1: Initialize ensemble and Hamiltonians"
 
    #----------------------- Initialization ---------------------------
    ens = Ensemble(ntraj, nel, nnucl) 

    ens.ham_set_ham("model", ham_indx)
    ens.ham_set_rep(rep)
    ens.ham_set_params([mass, omega, Er, V_, 1.5e-2+ie*dE ])


    THERM = Thermostat({"Q":100.0, "nu_therm":0.01, "NHC_size":2, "Temperature":300.0, "thermostat_type":"Nose-Hoover" })   
    THERM.set_Nf_t(1)
    THERM.set_Nf_r(0);
    THERM.init_nhc();
    therms = []
  
    for i in xrange(ntraj):
        # Array of thermostats
        therm = Thermostat({"Q":100.0, "nu_therm":0.01, "NHC_size":2, "Temperature":300.0, "thermostat_type":"Nose-Hoover" }) 
        therm.set_Nf_t(1);
        therm.set_Nf_r(0);
        therm.init_nhc();
        therms.append(therm)

        # Initialization of nuclear variables
        ens.mol[i].mass[0] = mass
        ens.mol[i].q[0], ens.mol[i].p[0] = boltz()
#        print i, ens.mol[i].q[0], ens.mol[i].p[0]
        ens.mol[i].f[0] = 0.0
        ens.ham_set_v(i, [ens.mol[i].p[0]/ens.mol[i].mass[0]] )
        ens.ham_set_q(i, ens.mol[i].q )


        # Initialization of electronic variables
        if rep==0:    # Diabatic - we start in left diabatic well 
            ens.el[i].istate = isurface
            ens.ham_compute(i)

        elif rep==1:  # Adiabatic
        # We start in left diabatic state - this corresponds to a mixture of two adiabatic states
        # The TD-SE populations are set from the transformation coefficients:
        # Later - set Hamiltonian parameters, if vary them
            ens.ham_compute(i)

            # Adiabatic energies
            E0 = ens.ham_H(i,0,0).real
            E1 = ens.ham_H(i,1,1).real
            # Extract diabatic energies:
            ens.ham_set_rep(i, 0)
            H0 = ens.ham_H(i,0,0).real 
            H1 = ens.ham_H(i,1,1).real
            V  = ens.ham_H(i,0,1).real
            ens.ham_set_rep(i, rep) # set representation to the original one

            #print E0, E1, H0, H1, V

            f0 = V/math.sqrt((H0 - E0)*(H0 - E0) + V*V)          # f0 - left diabat on 0
            g0 = (H0 - E0)/math.sqrt((H0 - E0)*(H0 - E0) + V*V)  # g0 - right diabat on 0
            f1 = V/math.sqrt((H0 - E1)*(H0 - E1) + V*V)          # f1 - left diabat on 1
            g1 = (H0 - E1)/math.sqrt((H0 - E1)*(H0 - E1) + V*V)  # g1 - right diabat on 1

            # cout<<"Probability matrix:\n";
            # cout<<" f0(left on 0)= "<<f0*f0<<"    f1(left on 1)= "<<f1*f1<<"\n";
            # cout<<" g0(right on 0)= "<<g0*g0<<"    g1(right on 1)= "<<g1*g1<<"\n";

            # Inverse C:
            # It is easy to prove that transpose is equal to inverse, so we don't need explicit inverse (below)
            # But still, lets use direct inverse:
            #den = g1*f0 - g0*f1;
            #a00 = g1/den;
            #a01 = -f1/den;
            #a10 = -g0/den;
            #a11 = f0/den;

            ksi = rnd.uniform(-2.0*math.pi, 2.0*math.pi);

            # We start on left diabat, so use only f (a01 and a11) coefficients, phases may be arbitrary
            ens.el[i].q[0] = f0 * math.cos(ksi)
            ens.el[i].p[0] = f0 * math.sin(ksi)

            ens.el[i].q[1] = f1 * math.cos(ksi)
            ens.el[i].p[1] = f1 * math.sin(ksi)

            # Set initial discrete state randomly:
            ksi = rnd.uniform(0.0, 1.0)
       
            if(ksi<f0*f0):
                ens.el[i].istate = 0  # 0-th adiabatic state
            else:
                ens.el[i].istate = 1
            #print i, ens.el[i], ens.el[i].istate

            


    # Update energies and forces for this setup
    epot = compute_potential_energy(ens, 1)  # 0 - MF forces, 1 - FSSH
    compute_forces(ens, 1)

    ekin = compute_kinetic_energy(ens)

    e_bath = 0.0
    if t_method=="nhc":
        e_bath = THERM.energy()
    elif t_method=="ens_nhc":
        for i in xrange(ntraj):
            e_bath = e_bath + therms[i].energy()
        e_bath = e_bath / float(ens.ntraj)


    e_tot = epot + ekin
    e_ext = e_tot + e_bath

    pops = ens.sh_pop1(-1000000, 1000000)
    x, sigma_ = x_ave(ens)
    snap = 0


    f = open("relax"+str(ie)+".txt","w")
    f.write("%10.8f   %10.8f   %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n" % (snap*nstep*dt, pops[0], pops[1], x, sigma_, e_tot, e_ext))
    f.close()


    #---------------------- Propagation -------------------------------    
    S = []
    for i in xrange(ens.ntraj):
        S.append(0.0)

    print "Step 2: Propagation"
    for snap in xrange(nsnap):
        for step in xrange(nstep):

            propagate_ensemble(dt, ens, 1)
        
            epot = compute_potential_energy(ens, 1)  # 0 - MF forces, 1 - FSSH
            ekin = compute_kinetic_energy(ens)
            e_tot = epot + ekin
        
        
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
                    old = ens.el[i].istate
                    ens.el[i].istate = hop(ens.el[i].istate, ens, i , ksi, g, do_rescaling, rep, do_reverse)  # this operation will also rescale velocities, if necessary

                    #print g.get(0,1)
                    #if(old!=ens.el[i].istate):
                    #    print snap, step, i
                    #print "%5i -> %5i" % (old, ens.el[i].istate)
        
        #---------------------- Compute and print info for given snap ------------------------------

        pops = ens.sh_pop1(-1000000, 1000000)
        x, sigma_ = x_ave(ens)

        e_bath = 0.0
        if t_method=="nhc":
            e_bath = THERM.energy()
        elif t_method=="ens_nhc":
            for i in xrange(ntraj):
                e_bath = e_bath + therms[i].energy()
            e_bath = e_bath / float(ens.ntraj)

        e_ext = e_tot + e_bath


        f = open("relax"+str(ie)+".txt","a")
        f.write("%10.8f   %10.8f   %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n" % ((snap+1)*nstep*dt, pops[0], pops[1], x, sigma_, e_tot, e_ext))
        f.close()

    


print "time propagation is done"      


