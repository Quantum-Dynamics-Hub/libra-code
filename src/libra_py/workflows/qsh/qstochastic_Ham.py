#*********************************************************************************
#* Copyright (C) 2018 Wei Li and Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
  Implementation of the Quasistochastic Hamiltonian method
      Akimov, J. Phys. Chem. Lett. 2017, 8, 5190
"""

import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
#from libra_py import *

import libra_py.workflows.common_utils as comn
import libra_py.units as units
import libra_py.acf_vector as acf_vector
import libra_py.tsh as tsh


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
        #print X[n].get(a,b)
    
    
    # Now compute ACFs of X matrix elements and print out the corresponding data
    T,  norm_acf,  raw_acf  = acf_vector.acf( acf_vector.center_data(data_ab)  , dt )  # dt is in a.u.
    
    dw =  1.0   # in cm^-1
    wspan = 3000.0 # in cm^-1
    dw = dw * units.inv_cm2Ha        # convert to Ha (atomic units)
    wspan = wspan * units.inv_cm2Ha        # convert to Ha (atomic units)
    
    f = open(filename+"_acf"+str(a)+"_"+str(b)+".txt","w")   
    tsz = len(T)
    for it in xrange(tsz):
      f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it]*units.au2fs , norm_acf[it], raw_acf[it]))
    f.close()
    
    # Do the FT
    W,  J  = acf_vector.ft(norm_acf,  wspan, dw, dt)  # dt is in a.u.
    jsz = len(W)
    
    f = open(filename+"_spectrum"+str(a)+"_"+str(b)+".txt","w")
    sp = MATRIX(jsz, 1)
    for iw in xrange(jsz):
        f.write("%8.5f  %8.5f  \n" % (W[iw]*units.au2wavn, J[iw] ) )
        sp.set(iw, J[iw]*J[iw])
    f.close()

    
    # Determine all frequencies (peaks) and sort them (in accending mannaer)
    out = comn.find_maxima(sp)

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
        freqs.append( [W[indx]*units.au2wavn, J[indx], J[indx]/norm ] )

        lgfile.write("index= %3i  frequency= %8.5f  amplitude= %8.5f normalized_amplitude= %8.5f \n" % (i, W[indx]*units.au2wavn, J[indx], J[indx]/norm) )
    lgfile.close()    



    
    lgfile = open("run.log", "a")
    print "max frequency for ", filename, " = ", freqs

    for a in freqs:
        print " Timescale is = ", a[0]/units.au2wavn, " Ha", " omega = E/hbar ", a[0]/units.au2wavn, " 2 pi*a.u. of time^-1",\
              " linear frequency = ", (a[0]/units.au2wavn)/(2.0*math.pi), " a.u.^-1", " Timescale = ", 2.0*math.pi*units.au2wavn/a[0], " a.u. of time ",\
              2.0*math.pi*units.au2wavn*units.au2fs/a[0], " fs, Amplitude = ", a[1], " Normalized amplitude = ", a[2]

#        lgfile.write("Timescale is = %8.5f Ha, omega = E/hbar %8.5f 2 pi*a.u. of time-1 linear frequence = %8.5f a.u-1,\
#                     Timescale = %8.5f a.u. of time \n" % ( a[0]/au2wavn, a[0]/au2wavn, a[0]/(au2wavn*2.0*math.pi), 2.0*math.pi*au2wavn/a[0],  )      )  

        lgfile.write("Timescale is = %8.5f Ha, omega = E/hbar %8.5f 2 pi*a.u. of time^-1 \
                      linear frequence = %8.5f a.u.^-1  Timescale = %8.5f a.u. of time \
                      %8.5f fs, Amplitude = %8.5f Normalized amplitude = %8.5f \n" % 
                     ( (a[0]/units.au2wavn), (a[0]/units.au2wavn), ((a[0]/units.au2wavn)/(2.0*math.pi)), (2.0*math.pi*units.au2wavn/a[0]), (2.0*math.pi*units.au2wavn*units.au2fs/a[0]), a[1], a[2] ) 
                    )

    lgfile.close()
    # Approximate the random data as:
    # f(t) = A * sin(omega * t + delta) +  B
   
    return freqs




def compute_Hvib(Nfreqs, freqs_re0, freqs_re1, freqs_im, t, 
                 H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                 H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                 dev):

    Hvib_stoch_re = MATRIX(2,2)
    Hvib_stoch_im = MATRIX(2,2)

    fu = [0.0, 0.0, 0.0]
    for i in xrange(Nfreqs):
        fu[0] = fu[0] + freqs_re0[i][2] * math.sin(freqs_re0[i][0]*t/units.au2wavn)
        fu[1] = fu[1] + freqs_re1[i][2] * math.sin(freqs_re1[i][0]*t/units.au2wavn)
        fu[2] = fu[2] + freqs_im[i][2] * math.sin(freqs_im[i][0]*t/units.au2wavn)
    
    xab = H_vib_re_ave.get(0,0) + H_vib_re_std.get(0,0) * (fu[0]/dev[0] ) 
    if xab < dw_Hvib_re.get(0,0):
        xab = dw_Hvib_re.get(0,0)
    elif xab > up_Hvib_re.get(0,0):
        xab = up_Hvib_re.get(0,0)
    Hvib_stoch_re.set(0,0,   xab )

    xab = H_vib_re_ave.get(1,1) + H_vib_re_std.get(1,1) * (fu[1]/dev[1] ) 
    if xab < dw_Hvib_re.get(1,1):
        xab = dw_Hvib_re.get(1,1)
    elif xab > up_Hvib_re.get(1,1):
        xab = up_Hvib_re.get(1,1)
    Hvib_stoch_re.set(1,1,   xab )

    xab = H_vib_im_ave.get(0,1) + H_vib_im_std.get(0,1) * (fu[2]/dev[2] ) 
    if xab < dw_Hvib_im.get(0,1):
        xab = dw_Hvib_im.get(0,1)
    elif xab > up_Hvib_im.get(0,1):
        xab = up_Hvib_im.get(0,1)
    Hvib_stoch_im.set(0,1,   xab )
    Hvib_stoch_im.set(1,0,  -xab )

    Hvib_stoch = CMATRIX(Hvib_stoch_re, Hvib_stoch_im)

    return Hvib_stoch



def run(params):

    use_boltz_factor = 1;
    dt = params["dt"]
    T = params["T"]
    norbitals = params["norbitals"]  # the number of orbitals in the input files
    act_sp = Py2Cpp_int(params["active_space"])
    nsteps = params["nsteps"]
    nfiles = params["nfiles"]
    nstates = len(act_sp)
    istate = params["istate"]
    rt = params["rt"]
    Nfreqs = params["Nfreqs"]
    ntraj = params["ntraj"]
    set_deco = params["set_deco"]
    deco_time = params["deco_time"] * units.fs2au

    rnd = Random()

    ################## Part 1: Collect files and compute statistics ================

    # Electronic Hamiltonian
    H_vib_re = []  # list of MATRIX
    H_vib_im = []  # list of MATRIX
    H_vib = [] # list of CMATRIX
    dH = [] # list of MATRIX
    dev = [0.0, 0.0, 0.0]

    H_vib_re_ave, H_vib_re_std = None, None
    H_vib_im_ave, H_vib_im_std = None, None
    dw_Hvib_re, up_Hvib_re = None, None
    dw_Hvib_im, up_Hvib_im = None, None

    freqs_re, freqs_im = None, None


    for i in xrange(0, nfiles): # how many files we have

        ##############################################################################
        # Read in the "elementary" overlaps and energies - in the basis of KS orbitals
        ##############################################################################       

        filename_re = rt+params["Hvib_re_prefix"]+str(i)+params["Hvib_re_suffix"]
        filename_im = rt+params["Hvib_im_prefix"]+str(i)+params["Hvib_im_suffix"]
        Hvib = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, act_sp)
        Hvib.scale(-1, -1, 0.5) #convert from Ry to Ha 
        H_vib.append(Hvib)

        hvib_re, hvib_im = Hvib.real(), Hvib.imag()
        H_vib_re.append(hvib_re)
        H_vib_im.append(hvib_im)

    
    H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re = comn.mat_stat(H_vib_re)
    H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im = comn.mat_stat(H_vib_im)    
    print "dw_re = "; dw_Hvib_re.show_matrix()
    print "up_re = "; up_Hvib_re.show_matrix()
    print "dw_im = "; dw_Hvib_im.show_matrix()
    print "up_im = "; up_Hvib_im.show_matrix()


    freqs_re0 = mat_freqs(H_vib_re, 0, 0, dt, "_H_vib_re_E0_", Nfreqs)
    freqs_re1 = mat_freqs(H_vib_re, 1, 1, dt, "_H_vib_re_E1_", Nfreqs)
    freqs_im = mat_freqs(H_vib_im, 0, 1, dt, "_H_vib_im_D01_", Nfreqs)

#    dH_ave, dH_std, dw_dH, up_dH =  comn.mat_stat(comn.energy_gaps(H_vib))

    # whether to set decoherence maunally
    if set_deco == 0:
       decoh_times, decoh_rates = comn.decoherence_times(Hvib, 1)

    elif set_deco ==1:
       decoh_times = MATRIX(nstates, nstates)
       decoh_rates = MATRIX(nstates, nstates)
       for a in xrange(nstates):
          for b in xrange(nstates):
             if a==b:
                decoh_times.set(a,a, 1000000.0)
                decoh_rates.set(a,a, 0.0)
             else:
                tau = deco_time # in a.u. unit
                decoh_times.set(a,b, tau)
                decoh_rates.set(a,b, 1.0/tau)

    if len(freqs_re0) < Nfreqs:
        Nfreqs = len(freqs_re0)
        print "The input Nfreqs is larger than the maximal number of the peaks, changing it to ", Nfreqs

    # Ok, now we have the function - sum of sines, so let's compute the standard deviation
    # This is a silly method - just do it numerically

    for r in xrange(1000000):

        fu = [0.0, 0.0, 0.0]
        for i in xrange(Nfreqs):
            fu[0] = fu[0] + freqs_re0[i][2] * math.sin(freqs_re0[i][0]*r*dt/units.au2wavn)
            fu[1] = fu[1] + freqs_re1[i][2] * math.sin(freqs_re1[i][0]*r*dt/units.au2wavn)
            fu[2] = fu[2] + freqs_im[i][2] * math.sin(freqs_im[i][0]*r*dt/units.au2wavn)

        for i in xrange(3):
            dev[i] = dev[i] + fu[i]**2

    for i in xrange(3):
        dev[i] = math.sqrt( dev[i] / 1000000.0 )
                    

    bf = math.exp( - ( H_vib_re_ave.get(1,1) - H_vib_re_ave.get(0,0) )/(units.kB*T) )
    print "Boltzmann factor = ", bf
    # bf = e/g,  e + g = 1;   e + e/bf = 1 =>  e *(1 + 1/bf) = 1 =>  e = 1/ (1 + 1/bf) = bf / (1 + bf)
    print "Equilibrium Ex st. population = ", bf /( 1.0 + bf )
    print freqs_re0, freqs_re1, freqs_im

    tau = [0.0, 0.0, 0.0]                 
    tau[0] = (2.0*math.pi*units.au2wavn/freqs_re0[0][0] )*dt
    tau[1] = (2.0*math.pi*units.au2wavn/freqs_re1[0][0] )*dt
    tau[2] = (2.0*math.pi*units.au2wavn/freqs_im[0][0] )*dt
 
    print tau, "in a.u."


    ################## Part 2: Run the dynamics with q-stochastic Hamiltonian ================


    out = open("populations.txt", "w")
    out.close()

 
    t = 0.0

    Cadi = []
    state = []

    for traj in xrange(ntraj):
        Cadi.append( CMATRIX(nstates,1) )
        Cadi[traj].set(istate, 0, 1.0+0.0j)
        state.append( istate  )
    

    Hvib = compute_Hvib(Nfreqs, freqs_re0, freqs_re1, freqs_im, t, 
                        H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                        H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                        dev)

    for i in xrange(1, nsteps): # nsteps     

        for traj in xrange(ntraj):
#            print traj
            propagate_electronic(0.5*dt, Cadi[traj], Hvib)
            Cadi[traj] = msdm(Cadi[traj], 0.5*dt, state[traj], decoh_rates)

    
        Hvib = compute_Hvib(Nfreqs, freqs_re0, freqs_re1, freqs_im, i*dt, 
                            H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                            H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                            dev)

        for traj in xrange(ntraj):
#            print traj
            propagate_electronic(0.5*dt, Cadi[traj], Hvib)
            Cadi[traj] = msdm(Cadi[traj], 0.5*dt, state[traj], decoh_rates)
    
        

        denmat_se = []
        denmat_sh = []
       
        for traj in xrange(ntraj):    
            g = compute_hopping_probabilities_fssh(Cadi[traj], Hvib, dt, use_boltz_factor, T);

            ksi = rnd.uniform(0.0,1.0);
            state[traj] = hop(state[traj], g, ksi);   # int
                  
            dm = Cadi[traj] * Cadi[traj].H() 

            denmat_se.append(dm)
            denmat_sh.append(dm)


        pop = tsh.update_sh_pop( state , nstates )
        dm_sh, dm_se = tsh.ave_pop(denmat_sh, denmat_se)


        out = open("populations.txt", "a")
        rec = (i*dt/41.0, dm_se.get(0,0).real, dm_se.get(1,1).real, pop[1], Hvib.get(0,0).real, Hvib.get(1,1).real, Hvib.get(0,1).imag  )
        out.write("%8.5f %8.5f %8.5f  %8.5f  %8.5f %8.5f %8.5f \n" %  rec )
        out.close()
      


        
