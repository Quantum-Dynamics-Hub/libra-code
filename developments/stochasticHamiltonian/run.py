#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

#####################################################################################
#
# This is the example code to run stochastic and quasi-stochastic Hamiltonian
# NA-MD calculations
#
#####################################################################################

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

au2fs = 0.02419 # 40 a.u. is 1 fs 
inv_cm2ev = (1.0/8065.54468111324)
ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
inv_cm2Ha = inv_cm2ev * ev2Ha
au2wavn = 27.211 * 8065.54468111324

def find_maxima(s):
# This function finds all the maxima of the data set and sorts them according to the data
# The maxima are defined as s[i-1] < s[i] > s[i+1]
# s - it the matrix of data 
#
# Returns the list of the indices of the maximal values

    max_indxs = []
    sz = s.num_of_elems
    for i in xrange(1, sz-1):
        if s.get(i) > s.get(i-1) and s.get(i) > s.get(i+1):
            max_indxs.append(i)

    inp = []
    sz = len(max_indxs)
    for i in xrange(sz):
        inp.append( [ max_indxs[i], s.get(max_indxs[i]) ] )


    out = merge_sort(inp)  # largest in the end


    lgfile = open("run.log", "a")
    lgfile.write("Found maxima of the spectrum:\n")
    for i in xrange(sz):
        lgfile.write("index = %3i  frequency index = %8.5f  intensity = %8.5f \n" % (i, out[sz-1-i][0], out[sz-1-i][1]) )
    lgfile.close()    
    
    return out
    



def flt_stat(X):
    N = len(X)
    res = 0.0

    #===== Average ====
    for i in xrange(N):
        res = res + X[i]
    res = res / float(N)


    #===== Std ========
    res2 = 0.0
    
    for i in xrange(N):
        res2 = res2 + (X[i] - res)**2
    res2 = math.sqrt( res2 / float(N) )

    return res, res2


def mat_stat(X):
# X - is a list of matrices

    N = len(X)
    res = MATRIX(X[0]); res *= 0.0

    #===== Average ====
    for i in xrange(N):
        res = res + X[i]
    res = res / float(N)


    #===== Std ========
    res2 = MATRIX(X[0]); res2 *= 0.0

    for a in xrange(res2.num_of_rows):
        for b in xrange(res2.num_of_cols):
        
            tmp = 0.0
            for i in xrange(N):
                tmp = tmp + (X[i].get(a,b) - res.get(a,b))**2
            tmp = math.sqrt( tmp / float(N) )

            res2.set(a,b, tmp)

    # Find maximal and minimal values
    up_bound = MATRIX(X[0]); up_bound *= 0.0
    dw_bound = MATRIX(X[0]); dw_bound *= 0.0



    for a in xrange(res2.num_of_rows):
        for b in xrange(res2.num_of_cols):

            up_bound.set(a,b, X[0].get(a,b))
            dw_bound.set(a,b, X[0].get(a,b))

            for i in xrange(N):
                xab = X[i].get(a,b)
                if xab > up_bound.get(a,b):
                    up_bound.set(a,b, xab)
                if xab < dw_bound.get(a,b):
                    dw_bound.set(a,b,xab)


    return res, res2, dw_bound, up_bound



def mat_freqs(X, a, b, dt, filename, Nfreqs):
# X - is a list of matrices
# a, b - indices that define which matrix elements to analyze
# dt - time step in a.u.
# filename - prefix for the filename to which the data will be printed out
# Nfreqs - the number of frequencies we want to extract

    N = len(X)
    sz = X[0].num_of_rows
    freqs = []
            

    # Collect info in a different format
    data_ab = []
    for n in xrange(N):
        data_ab.append(VECTOR(X[n].get(a,b), 0.0, 0.0))
    
    # Now compute ACFs of X matrix elements and print out the corresponding data
    T,  norm_acf,  raw_acf  = acf.acf( acf.center_data(data_ab)  , dt )  # dt is in a.u.
    
    dw =  1.0   # in cm^-1
    wspan = 3000.0 # in cm^-1
    dw = dw * inv_cm2Ha        # convert to Ha (atomic units)
    wspan = wspan * inv_cm2Ha        # convert to Ha (atomic units)
    
    f = open(filename+"_acf"+str(a)+"_"+str(b)+".txt","w")   
    tsz = len(T)
    for it in xrange(tsz):
      f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it]*au2fs , norm_acf[it], raw_acf[it]))
    f.close()
    
    # Do the FT
    W,  J  = acf.ft(norm_acf,  wspan, dw, dt)  # dt is in a.u.
    jsz = len(W)
    
    f = open(filename+"_spectrum"+str(a)+"_"+str(b)+".txt","w")
    sp = MATRIX(jsz, 1)
    for iw in xrange(jsz):
        f.write("%8.5f  %8.5f  \n" % (W[iw]*au2wavn, J[iw] ) )
        sp.set(iw, J[iw]*J[iw])
    f.close()

    
    # Determine all frequencies (peaks) and sort them (in accending mannaer)
    out = find_maxima(sp)

    #print out

    lgfile = open("run.log", "a")
    lgfile.write("Maximal peaks in the file "+filename+"_spectrum"+str(a)+"_"+str(b)+".txt\n")

    if Nfreqs > len(out):
        Nfreqs = len(out)

    szo = len(out) - 1

    norm = 0.0
    for i in xrange(Nfreqs):
        indx = out[szo-i][0]
        norm = norm + abs(J[indx])
       
    for i in xrange(Nfreqs):
        indx = out[szo-i][0]
        freqs.append( [W[indx]*au2wavn, J[indx], J[indx]/norm ] )

        lgfile.write("index= %3i  frequency= %8.5f  amplitude= %8.5f normalized_amplitude= %8.5f \n" % (i, W[indx]*au2wavn, J[indx], J[indx]/norm) )
    lgfile.close()    



    
    lgfile = open("run.log", "a")
    print "max frequency for ", filename, " = ", freqs

    for a in freqs:
        print " Timescale is = ", a[0]/au2wavn, " Ha", " omega = E/hbar ", a[0]/au2wavn, " 2 pi*a.u. of time^-1",\
              " linear frequency = ", (a[0]/au2wavn)/(2.0*math.pi), " a.u.^-1", " Timescale = ", 2.0*math.pi*au2wavn/a[0], " a.u. of time ",\
              2.0*math.pi*au2wavn*au2fs/a[0], " fs, Amplitude = ", a[1], " Normalized amplitude = ", a[2]

#        lgfile.write("Timescale is = %8.5f Ha, omega = E/hbar %8.5f 2 pi*a.u. of time-1 linear frequence = %8.5f a.u-1,\
#                     Timescale = %8.5f a.u. of time \n" % ( a[0]/au2wavn, a[0]/au2wavn, a[0]/(au2wavn*2.0*math.pi), 2.0*math.pi*au2wavn/a[0],  )      )  

        lgfile.write("Timescale is = %8.5f Ha, omega = E/hbar %8.5f 2 pi*a.u. of time^-1 \
                      linear frequence = %8.5f a.u.^-1  Timescale = %8.5f a.u. of time \
                      %8.5f fs, Amplitude = %8.5f Normalized amplitude = %8.5f \n" % 
                     ( (a[0]/au2wavn), (a[0]/au2wavn), ((a[0]/au2wavn)/(2.0*math.pi)), (2.0*math.pi*au2wavn/a[0]), (2.0*math.pi*au2wavn*au2fs/a[0]), a[1], a[2] ) 
                    )

    lgfile.close()
    # Approximate the random data as:
    # f(t) = A * sin(omega * t + delta) +  B
    # 
    #
   
    return freqs




rnd = Random()
kb = 3.166811429e-6

##============== Simulation parameters ================
nsnap = 40000
nstep = 41
dt = 1.0  # a.u.

therm_size = int(1.1 * nsnap * nstep)   # for how long to run NVT without SH
sampl_size = int(0.05 * nsnap * nstep ) # for how long to run explicit Hamiltonian calculations

Ntraj = 1000
use_boltz_factor = 1
do_rescaling = 0
do_reverse = 0


# A flag to control the construction of the random Hamiltonian
# 0 : mean + random fluctuations
# 1 : mean + amplitude * sin(a number of modes)
rand_Ham = 1

# Correlate random numbers in stochastic Ham generation:
# 0 - do not correlate
# 1 - do correlate
corr_opt = 0

# Randon fluctuations in stochastic Ham. generation:
# 0 - uniformely distributed
# 1 - normal distribution
ksi_opt = 1

# How many frequencies we want to use in the computations of the
# deterministic stochastic Hamiltonian
# if this number is larger than the actual number of the found frequencies,
# the latter will be used.
Nfreqs = 200

# Option on if we want to restrict the random (or deterministic) values of the
# energies (but not NACs) to a given window of the values - derived from the
# sampling
# 0 - do not restrict, use values as they are
# 1 - do restrict
do_filtering = 1




nu_therm = 0.001
T = 300.0

ham_indx = 7
E0 = -1.0 * kb * T
E1 = 1.0 * kb * T

V01=  1.0 * kb * T  #0.001 # * math.sqrt(10.0)
D = 1.0 * kb * T  # 0.001
print "E0 = ", E0, " E1 = ", E1, " D = ", D, " V01 = ", V01
L = 5.0
a = 0.1

lgfile = open("run.log", "w")
lgfile.write("Simulation parameters\n")
lgfile.write("nsnap = %6i \n" % nsnap)
lgfile.write("nstep = %6i \n" % nstep)
lgfile.write("dt (a.u.) = %8.5f \n" % dt)
lgfile.write("nsnap * nstep = %8i \n" % nsnap * nstep)
lgfile.write("therm_size (steps) = %6i \n" % therm_size)
lgfile.write("sampl_size (steps) = %6i \n" % sampl_size)
lgfile.write("Ntraj = %5i \n" % Ntraj)
lgfile.write("use_boltz_factor = %3i \n" % use_boltz_factor)
lgfile.write("do_rescaling = %3i \n" % do_rescaling)
lgfile.write("do_reverse = %3i \n" % do_reverse)
lgfile.write("corr_opt = %3i \n" % corr_opt)
lgfile.write("ksi_opt = %3i \n" % ksi_opt)
lgfile.write("rand_Ham = %3i \n" % rand_Ham)
lgfile.write("nu_therm = %8.5f \n" % nu_therm)
lgfile.write("T = %8.5f \n" % T)
lgfile.write("ham_indx = %3i \n" % ham_indx)
lgfile.write("E0 = %8.5f  E0(in kT) = %8.5f \n" % (E0, E0/(kb*T) ) )
lgfile.write("E1 = %8.5f  E1(in kT) = %8.5f \n" % (E1, E1/(kb*T) ) )
lgfile.write("V01 = %8.5f  V01(in kT) = %8.5f \n" % (V01, V01/(kb*T) ) )
lgfile.write("D = %8.5f  D(in kT) = %8.5f \n" % (D, D/(kb*T) ) )
lgfile.write("L = %8.5f  \n" % L )
lgfile.write("a = %8.5f  \n" % a )
lgfile.close()

##=======================================================


# First, create the Hamiltonian
ham = Hamiltonian_Model(ham_indx)  # periodic SIN potential
ham.set_rep(1)  # adiabatic
ham.set_params([E0, E1, V01, D, L ])


# Here goes the external Hamiltonian - the one which will actually be used
ham_ex = Hamiltonian_Extern(2,1)  # 2 - # of electronic DOF, 1 - # of nuclear DOF
ham_ex.set_rep(1)  # adiabatic
ham_ex.set_adiabatic_opt(0)  # use the externally-computed adiabatic electronic Hamiltonian and derivatives
ham_ex.set_vibronic_opt(0)  # use the externally-computed vibronic Hamiltonian and derivatives

# Actual matrices
ham_adi = MATRIX(2,2);  ham_ex.bind_ham_adi(ham_adi);
d1ham_adi = MATRIXList()

tmp = MATRIX(2,2)
d1ham_adi.append(tmp);  ham_ex.bind_d1ham_adi(d1ham_adi);

ham_vib = CMATRIX(2,2);  ham_ex.bind_ham_vib(ham_vib);



THERM = Thermostat({"nu_therm":nu_therm, "NHC_size":10, "Temperature":T, "thermostat_type":"Nose-Hoover" })   
THERM.set_Nf_t(1)
THERM.set_Nf_r(0);
THERM.init_nhc();



# Electronic - 2 levels, starting at 1-th state
el = Electronic(2,1)
el.istate = 0         # originally start at the GS (to run thermalization dynamics and sample Hamiltonian according to NBRA)

# Nuclear DOFs
mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = -10.0 
mol.p[0] = math.sqrt( kb * T * mol.mass[0] ) 

lgfile = open("run.log", "a")
lgfile.write("mass = %8.5f \n" % mol.mass[0])
lgfile.write("q[0] = %8.5f \n" % mol.q[0])
lgfile.write("p[0] = %8.5f \n" % mol.p[0])
lgfile.close()


f = open("_tsh_extern.txt","w")
f.close()


# Update matrices
ham.set_v([ mol.p[0]/mol.mass[0] ])
ham.set_q([ mol.q[0] ])
ham.compute()
for i in [0,1]:
    for j in [0,1]:
        ham_adi.set(i,j, ham.H(i,j).real)
        ham_vib.set(i,j, ham.Hvib(i,j))
        d1ham_adi[0].set(i,j, ham.dHdq(i,j,0).real)


# Initialization
ham_ex.set_v([ mol.p[0]/mol.mass[0] ])
epot = compute_forces(mol, el, ham_ex, 1)   # FSSH forces
ekin = compute_kinetic_energy(mol) 
ebath = THERM.energy()


sh_pop = [0.0, 1.0]
states = [1]*Ntraj    # but initiate all trajectories to start at the highest energy state


##==========================================
# Stochastic variables:

# Electronic Hamiltonian
H_vib_re = []  # list of MATRIX
H_vib_im = []  # list of MATRIX
p0 = [] # list of double

H_vib_re_ave, H_vib_re_std = None, None
H_vib_im_ave, H_vib_im_std = None, None
dw_Hvib_re, up_Hvib_re = None, None
dw_Hvib_im, up_Hvib_im = None, None

p0_ave, p0_std = None, None

freqs_re, freqs_im = None, None
ksi = [0.0, 0.0, 0.0]
zeta = [0.0, 0.0, 0.0]
f_correl = [0.0, 0.0, 0.0]
f_correl2 = [0.0, 0.0, 0.0]
ksi_prev = [0.0, 0.0, 0.0]

##==========================================

count = 0  # counter of the steps (step size dt)

dev = [0.0, 0.0, 0.0]

for snap in xrange(nsnap):

    t = snap*nstep*dt

    f = open("_tsh_extern.txt","a")
    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  \
          etot= %8.5f  eext= %8.5f |c0|^2= %8.5f  |c1|^2= %8.5f  Re|c01|= %8.5f istate= %8.5f \
          sh_pop0= %8.5f sh_pop1= %8.5f, E0= %8.5f E1= %8.5f d= %8.5f, \
          E0= %8.5f E1= %8.5f d= %8.5f \n" % 
            (t, mol.q[0], mol.p[0], ekin, epot, 
          ekin+epot, ekin+epot+ebath, el.rho(0,0).real, el.rho(1,1).real, el.rho(0,1).real, el.istate, sh_pop[0], sh_pop[1],
             ham_ex.Hvib(0,0).real, ham_ex.Hvib(1,1).real, ham_ex.Hvib(0,1).imag,
             ham_vib.get(0,0).real, ham_vib.get(1,1).real, ham_vib.get(0,1).imag)
#            ) 
           )
    f.close()


    for step in xrange(nstep):

        t = (snap*nstep + step)*dt

        #================= Propagation of electrons and nuclei ===========================
        # el-dyn
        el.propagate_electronic(0.5*dt, ham_ex)
        
 
        if count <= sampl_size:  # explicit Hamiltonian calculations

            #============== Regular v-Verlet =============
            # Nuclear dynamics
            if count<=therm_size:
                mol.p[0] = mol.p[0] * THERM.vel_scale(0.5*dt)

            mol.propagate_p(0.5*dt)        
            mol.propagate_q(dt)
        
            # Update matrices
            ham.set_v([ mol.p[0]/mol.mass[0] ])
            ham.set_q([ mol.q[0] ])
            ham.compute()

            for i in [0,1]:
                for j in [0,1]:
                    ham_adi.set(i,j, ham.H(i,j).real)
                    ham_vib.set(i,j, ham.Hvib(i,j))
                    d1ham_adi[0].set(i,j, ham.dHdq(i,j,0).real)


            epot = compute_forces(mol, el, ham_ex, 1)  # FSSH forces
            ekin = compute_kinetic_energy(mol) 

            if count<=therm_size:        
                THERM.propagate_nhc(dt, ekin, 0.0, 0.0)
        
            mol.propagate_p(0.5*dt)

            if count<=therm_size:
                mol.p[0] = mol.p[0] * THERM.vel_scale(0.5*dt)

            ekin = compute_kinetic_energy(mol)
            ebath = THERM.energy()

            #==================================================
            #===== Sampling =========
            hel = MATRIX(2,2)
            hvib= CMATRIX(2,2)
            for i in [0,1]:
                for j in [0,1]:
                    hel.set(i,j, ham.H(i,j).real)
                    hvib.set(i,j, ham.Hvib(i,j))
            H_vib_re.append(hvib.real())
            H_vib_im.append(hvib.imag())
            p0.append(mol.p[0])

            #======= Generate distributions =======
            if count==sampl_size:

                H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re = mat_stat(H_vib_re)
                H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im = mat_stat(H_vib_im)

          
                print "dw_re = "; dw_Hvib_re.show_matrix()
                print "up_re = "; up_Hvib_re.show_matrix()
                print "dw_im = "; dw_Hvib_im.show_matrix()
                print "up_im = "; up_Hvib_im.show_matrix()


                lgfile = open("run.log", "a")
                lgfile.write("dw_re = \n")
                sz1 = dw_Hvib_re.num_of_rows
                sz2 = dw_Hvib_re.num_of_cols
                for a in xrange(sz1):
                    line = ""
                    for b in xrange(sz2):
                        line = line + " %8.5f " % dw_Hvib_re.get(a,b)
                    line = line + "\n"
                    lgfile.write(line)
                lgfile.write("\n")
                lgfile.close()


                lgfile = open("run.log", "a")
                lgfile.write("up_re = \n")
                for a in xrange(sz1):
                    line = ""
                    for b in xrange(sz2):
                        line = line + " %8.5f " % up_Hvib_re.get(a,b)
                    line = line + "\n"
                    lgfile.write(line)
                lgfile.write("\n")
                lgfile.close()

                lgfile = open("run.log", "a")
                lgfile.write("dw_im = \n")
                for a in xrange(sz1):
                    line = ""
                    for b in xrange(sz2):
                        line = line + " %8.5f " % dw_Hvib_im.get(a,b)
                    line = line + "\n"
                    lgfile.write(line)
                lgfile.write("\n")
                lgfile.close()

                lgfile = open("run.log", "a")
                lgfile.write("up_im = \n")
                for a in xrange(sz1):
                    line = ""
                    for b in xrange(sz2):
                        line = line + " %8.5f " % up_Hvib_im.get(a,b)
                    line = line + "\n"
                    lgfile.write(line)
                lgfile.write("\n")
                lgfile.close()



                if rand_Ham==1 or (rand_Ham==0 and corr_opt==1):
                    freqs_re0 = mat_freqs(H_vib_re, 0, 0, dt, "H_vib_re_E0_", Nfreqs)
                    freqs_re1 = mat_freqs(H_vib_re, 1, 1, dt, "H_vib_re_E1_", Nfreqs)
                    freqs_im = mat_freqs(H_vib_im, 0, 1, dt, "H_vib_im_D01_", Nfreqs)


                    if len(freqs_re0) < Nfreqs:
                        Nfreqs = len(freqs_re0)
                        print "The input Nfreqs is larger than the maximal number of the peaks, changing it to ", Nfreqs


                # Ok, now we have the function - sum of sines, so let's compute the standard deviation
                # This is a silly method - just do it numerically

                for r in xrange(1000000):

                    fu = [0.0, 0.0, 0.0]
                    for i in xrange(Nfreqs):
                        fu[0] = fu[0] + freqs_re0[i][2] * math.sin(freqs_re0[i][0]*r*dt/au2wavn)
                        fu[1] = fu[1] + freqs_re1[i][2] * math.sin(freqs_re1[i][0]*r*dt/au2wavn)
                        fu[2] = fu[2] + freqs_im[i][2] * math.sin(freqs_im[i][0]*r*dt/au2wavn)

                    for i in xrange(3):
                        dev[i] = dev[i] + fu[i]**2

                for i in xrange(3):
                    dev[i] = math.sqrt( dev[i] / 1000000.0 )
                    
                  


                bf = math.exp( - ( H_vib_re_ave.get(1,1) - H_vib_re_ave.get(0,0) )/(kb*T) )
                print "Boltzmann factor = ", bf
                # bf = e/g,  e + g = 1;   e + e/bf = 1 =>  e *(1 + 1/bf) = 1 =>  e = 1/ (1 + 1/bf) = bf / (1 + bf)
                print "Equilibrium Ex st. population = ", bf /( 1.0 + bf )

                #freqs_re = [314.00000000000006, 316.0]
                #freqs_im = [333.00000000000006]

                print freqs_re0, freqs_re1, freqs_im
                p0_ave, p0_std = flt_stat(p0)

                mol.p[0] = 0.0
                mol.f[0] = 0.0

                lgfile = open("run.log", "a")
                lgfile.write("Boltzmann factor = %8.5f \n" % bf)
                lgfile.write("Equilibrium Ex st. population = %8.5f \n" %  (bf/( 1.0 + bf )) )

                
                lgfile.write("freqs_re_E0_: ")
                for it in freqs_re0:
                    lgfile.write(" %8.5f %8.5f %8.5f" %  (it[0], it[1], it[2]) )
                lgfile.write("\n")

                lgfile.write("freqs_re_E1_: ")
                for it in freqs_re1:
                    lgfile.write(" %8.5f %8.5f %8.5f" %  (it[0], it[1], it[2]) )
                lgfile.write("\n")


                lgfile.write("freqs_im_D01_: ")
                for it in freqs_im:
                    lgfile.write(" %8.5f %8.5f %8.5f" %  (it[0], it[1], it[2]) )
                lgfile.write("\n")

                lgfile.write("p0_ave = %8.5f p0_std = %8.5f \n" % (p0_ave, p0_std) )

                lgfile.close()

        
 
        else:  # Use the stochastic Hamiltonian 
 
              
            tau = [0.0, 0.0, 0.0]                 
            tau[0] = (2.0*math.pi*au2wavn/freqs_re0[0][0] )*dt
            tau[1] = (2.0*math.pi*au2wavn/freqs_re1[0][0] )*dt
            tau[2] = (2.0*math.pi*au2wavn/freqs_im[0][0] )*dt

            if count == sampl_size+1:
                lgfile = open("run.log", "a")
                lgfile.write("tau[0] = %8.5f tau[1] = %8.5f tau[2] = %8.5f\n" % (tau[0], tau[1], tau[2]))
                lgfile.close()

            for ind in [0,1,2]:
                if corr_opt==0:
                    f_correl[ind] = 0.0
                    f_correl2[ind] = 1.0
                elif corr_opt==1:
                    f_correl[ind] = math.exp(-1.0/tau[ind])
                    f_correl2[ind] = math.sqrt(1.0 - f_correl[ind]**2)

            if count == sampl_size+1:
                lgfile = open("run.log", "a")
                lgfile.write("f_correl[0] = %8.5f f_correl[1] = %8.5f f_correl[2] = %8.5f\n" % (f_correl[0], f_correl[1], f_correl[2]) )
                lgfile.write("f_correl2[0] = %8.5f f_correl2[1] = %8.5f f_correl2[2] = %8.5f\n" % (f_correl2[0], f_correl2[1], f_correl2[2]) )
                lgfile.close()


            if ksi_opt==0:
                ksi = [rnd.uniform(-1.0,1.0), rnd.uniform(-1.0,1.0), rnd.uniform(-1.0,1.0)]
            elif ksi_opt==1:
                ksi = [rnd.normal(), rnd.normal(), rnd.normal()]

            # Correlation part:
            zeta = [0.0, 0.0, 0.0]
            if count == sampl_size + 1:
                for ind in [0,1,2]:
                    zeta[ind] = ksi[ind]
            else:
                for ind in [0,1,2]:
                    zeta[ind] = f_correl[ind] * ksi_prev[ind] + f_correl2[ind] * ksi[ind]



            Hvib_stoch_re = MATRIX(2,2)
            Hvib_stoch_im = MATRIX(2,2)

            if rand_Ham==0:
                # Use random matrixes, no longer need forces.            
                #mol.p[0] = p0_ave + p0_std * zeta0
                #ekin = compute_kinetic_energy(mol) 

                xab = H_vib_re_ave.get(0,0) + H_vib_re_std.get(0,0) * zeta[0]
                if do_filtering==1:
                    if xab < dw_Hvib_re.get(0,0):
                        xab = dw_Hvib_re.get(0,0)
                    elif xab > up_Hvib_re.get(0,0):
                        xab = up_Hvib_re.get(0,0)

                Hvib_stoch_re.set(0,0,   xab )


                xab = H_vib_re_ave.get(1,1) + H_vib_re_std.get(1,1) * zeta[1]
                if do_filtering==1:
                    if xab < dw_Hvib_re.get(1,1):
                        xab = dw_Hvib_re.get(1,1)
                    elif xab > up_Hvib_re.get(1,1):
                        xab = up_Hvib_re.get(1,1)

                Hvib_stoch_re.set(1,1,  xab )

                xab = H_vib_im_ave.get(0,1) + H_vib_im_std.get(0,1) * zeta[2]
                Hvib_stoch_im.set(0,1,  xab )
                Hvib_stoch_im.set(1,0, -xab )


            elif rand_Ham==1:

#                mol.p[0] = p0_ave + p0_std * zeta[0]
#                THERM.propagate_nhc(dt, ekin, 0.0, 0.0)
#                ebath = THERM.energy()

                
                fu = [0.0, 0.0, 0.0]
                for i in xrange(Nfreqs):
                    fu[0] = fu[0] + freqs_re0[i][2] * math.sin(freqs_re0[i][0]*t/au2wavn)
                    fu[1] = fu[1] + freqs_re1[i][2] * math.sin(freqs_re1[i][0]*t/au2wavn)
                    fu[2] = fu[2] + freqs_im[i][2] * math.sin(freqs_im[i][0]*t/au2wavn)
      

                xab = H_vib_re_ave.get(0,0) + H_vib_re_std.get(0,0) * (fu[0]/dev[0] ) # + 0.1*zeta[0] )

                if do_filtering==1:
                    if xab < dw_Hvib_re.get(0,0):
                        xab = dw_Hvib_re.get(0,0)
                    elif xab > up_Hvib_re.get(0,0):
                        xab = up_Hvib_re.get(0,0)
                Hvib_stoch_re.set(0,0,   xab )


                xab = H_vib_re_ave.get(1,1) + H_vib_re_std.get(1,1) * (fu[1]/dev[1] ) # + 0.1*zeta[1] )

                if do_filtering==1:
                    if xab < dw_Hvib_re.get(1,1):
                        xab = dw_Hvib_re.get(1,1)
                    elif xab > up_Hvib_re.get(1,1):
                        xab = up_Hvib_re.get(1,1)
                Hvib_stoch_re.set(1,1,   xab )


                xab = H_vib_im_ave.get(1,1) + H_vib_im_std.get(0,1) * (fu[2]/dev[2] ) # + 0.1*zeta[2] )
#                if xab < dw_Hvib_im.get(0,1):
#                    xab = dw_Hvib_im.get(0,1)
#                elif xab > up_Hvib_im.get(0,1):
#                    xab = up_Hvib_im.get(0,1)
                Hvib_stoch_im.set(0,1,   xab )
                Hvib_stoch_im.set(1,0,  -xab )



                
            Hvib_stoch = CMATRIX(Hvib_stoch_re, Hvib_stoch_im)
            for i in [0,1]:
                for j in [0,1]:
                    ham_adi.set(i,j, Hvib_stoch_re.get(i,j))
                    ham_vib.set(i,j, Hvib_stoch.get(i,j))

            ham_ex.compute()

            #epot = compute_forces(mol, el, ham_ex, 1)  # FSSH forces
            #ekin = compute_kinetic_energy(mol) 


            # Current becomes new
            for ind in [0,1,2]:
                ksi_prev[ind] = ksi[ind]


        
        # el-dyn
        el.propagate_electronic(0.5*dt, ham_ex)
        
        
                
        
        #================= Now, incorporate surface hop ===========================
        g = MATRIX(2,2)
        
        #if count<=therm_size:
        
        # Just choose the TSH scheme below
        #compute_hopping_probabilities_mssh(mol, el, ham, g, dt, use_boltz_factor, T)
        compute_hopping_probabilities_fssh(mol, el, ham_ex, g, dt, use_boltz_factor, T)
        #compute_hopping_probabilities_gfsh(mol, el, ham, g, dt, use_boltz_factor, T)
        
        for traj in xrange(Ntraj):
            ksi = rnd.uniform(0.0, 1.0)
            states[traj] = hop(states[traj], g, ksi)  
            #states[traj] = hop(states[traj], mol, ham_ex, ksi, g, do_rescaling, 1, do_reverse)  # this operation will also rescale velocities, if necessary
            el.istate = 0

        #============================================================================

        count = count + 1


    sh_pop = tsh.update_sh_pop(states, 2)
