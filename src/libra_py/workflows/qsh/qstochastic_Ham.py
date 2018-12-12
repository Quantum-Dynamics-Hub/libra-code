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

import libra_py.workflows.common_utils as comn
import libra_py.units as units
import libra_py.acf_vector as acf_vector
import libra_py.tsh as tsh


def decoh_method(set_decoherence,nstates,Hvib,deco_time):
    """
    set_decoherence - possible value, -1, 0, 1
    nstates - number of states
    Hvib, nstates x nstates matrix, vibronic Ham
    deco_time - float, decoherence time in case we want to set it muanlly

    Returns: nstates x nstates matrix, decoherence time, decoherence rates
    """
    if set_decoherence == 0:
        # set decoherence time from the Libra module
        # Akimov and Prezhdo, J. Phys. Chem. Lett., 2013, 4, 3857
        decoh_times, decoh_rates = comn.decoherence_times(Hvib, 1)

    else:
       decoh_times = MATRIX(nstates, nstates)
       decoh_rates = MATRIX(nstates, nstates)

       for a in xrange(nstates):
           for b in xrange(nstates):
               if a==b:
                   decoh_times.set(a,a, 1000000.0)
                   decoh_rates.set(a,a, 0.0)
               else:
                   if set_decoherence ==1:
                       # set the decoherence time manually
                       tau = deco_time # in a.u. unit
                       decoh_times.set(a,b, tau)
                       decoh_rates.set(a,b, 1.0/tau)

                   elif set_decoherence == -1:
                       # don't include decoherence effect, set it to infinite
                       decoh_times.set(a,b, 1000000.0)
                       decoh_rates.set(a,b, 0.0)

    return decoh_times, decoh_rates

def mat_freqs(X, a, b, dt, filename, Nfreqs):
    """
    X - is a list of matrices
    a, b - indices that define which matrix elements to analyze
    dt - time step in a.u.
    filename - prefix for the filename to which the data will be printed out
    Nfreqs - the number of frequencies we want to extract

    Return: list of list, frequencies
    """
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
    
    f = open(filename+"_acf_"+str(a)+"_"+str(b)+".txt","w")   
    tsz = len(T)
    for it in xrange(tsz):
      f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it]*units.au2fs , norm_acf[it], raw_acf[it]))
    f.close()
    
    # Do the FT
    W,  J  = acf_vector.ft(norm_acf,  wspan, dw, dt)  # dt is in a.u.
    jsz = len(W)
    
    f = open(filename+"_spectrum_"+str(a)+"_"+str(b)+".txt","w")
    sp = MATRIX(jsz, 1)
    for iw in xrange(jsz):
        f.write("%8.5f  %8.5f  \n" % (W[iw]*units.au2wavn, J[iw] ) )
        sp.set(iw, J[iw]*J[iw])
    f.close()

    
    # Determine all frequencies (peaks) and sort them (in accending manner)
    out = comn.find_maxima(sp)

    #print out

    lgfile = open("run.log", "a")
    lgfile.write("Maximal peaks in the file "+filename+"_spectrum_"+str(a)+"_"+str(b)+".txt\n")

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




def compute_Hvib(Nfreqs, freqs, t, nstates, 
                 H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                 H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                 dev):
    """
    Compute the QSH Hamiltonians
    Nfreqs - the number of frequencies we want to extract
    freqs - 2D list contains the frequencies for energies and couplings
    t - time step in a.u.
    nstates - number of states
    H_vib_re_ave - nstates x nstates matrix, average of energies for direct Hamiltonian
    H_vib_im_ave - nstates x nstates matrix, average of couplings for direct Hamiltonian
    H_vib_re_std - nstates x nstates matrix, std of energies for direct Hamiltonian
    H_vib_im_std - nstates x nstates matrix, std of couplings for direct Hamiltonian
    up_Hvib_re - nstates x nstates matrix, maximum value of energies for direct Hamiltonian
    up_Hvib_im - nstates x nstates matrix, maximum value of couplings for direct Hamiltonian
    dw_Hvib_re - nstates x nstates matrix, minimal value of energies for direct Hamiltonian
    dw_Hvib_im - nstates x nstates matrix, minimal value of couplings for direct Hamiltonian
    std - nstates x nstates matrix, std for direct Hamiltonian
    
    Return: nstates x nstates complex matrix contains QSH Hamiltonian and couplings 
    """
    Hvib_stoch_re = MATRIX(nstates,nstates)
    Hvib_stoch_im = MATRIX(nstates,nstates)

    fu = [ [ 0.0 for i in xrange(nstates)] for j in xrange(nstates)]

    for i in xrange(nstates):
        for j in xrange(nstates):
            for k in xrange(Nfreqs):
		    fu[i][j] = fu[i][j] + freqs[i][j][k][2] * math.sin(freqs[i][j][k][0]*t/units.au2wavn)

    for i in xrange(nstates):
        for j in xrange(nstates):
	    if i==j:
	        xab = H_vib_re_ave.get(i,j) + H_vib_re_std.get(i,j) * (fu[i][j]/dev[i][j] ) 
                if xab < dw_Hvib_re.get(i,j):
                    xab = dw_Hvib_re.get(i,j)
                elif xab > up_Hvib_re.get(i,j):
                    xab = up_Hvib_re.get(i,j)
                Hvib_stoch_re.set(i,j,   xab )

            elif i<j:
		xab = H_vib_im_ave.get(i,j) + H_vib_im_std.get(i,j) * (fu[i][j]/dev[i][j] ) 
                if xab < dw_Hvib_im.get(i,j):
                    xab = dw_Hvib_im.get(i,j)
                elif xab > up_Hvib_im.get(i,j):
                    xab = up_Hvib_im.get(i,j)
                Hvib_stoch_im.set(i,j,   xab )
                Hvib_stoch_im.set(j,i,  -xab )

    Hvib_stoch = CMATRIX(Hvib_stoch_re, Hvib_stoch_im)

    return Hvib_stoch



def run_namd(params):
    """
    The main procedure to run NA-MD calculations according to QSH workflow
    
    === Required parameter keys:  ===

    params["nsteps"]           [int] - define the length of QSH-NAMD simulation
    params["active_space"]     [list of ints] - which orbitals we care about (indexing starts with 0)
    params["norbitals"]        [int] - how many lines/columns in the file
    params["nfiles"]           [int] - how many files to read
    params["qsh_Ham_prefix"]   [string] - the prefix of output files
    params["time_inteval"]     [int] - define the time interval between adjacent output files  
    params["deco_time"]        [float] - decoherenc time in case we want to set it manually

    params["dt"]               [double] - nuclear dynamics integration time step [in a.u. of time, default: 41.0]
    params["istate"]           [int] - index of the initial state [default: 0]  
    params["ntraj"]            [int] - the number of stochastic surface hopping trajectories [default: 1]
    params["T"]                [double] - temperature of nuclear/electronic dynamics [in K, default: 300.0]
    params["set_decoherence"]  [int] - selection of how we include decoherence effect

    """

    critical_params = ["norbitals", "active_space", "nfiles",  "nsteps","qsh_Ham_prefix" ]
    default_params = { "T":300.0, "ntraj":100, "Nfreqs":1, 
                       "set_decoherence":-1, "deco_time":100000, "dt":41.0, 
                       "istate":1, "time_inteval":1 }
    comn.check_input(params, default_params, critical_params)

    use_boltz_factor = 1;
    dt = params["dt"]
    T = params["T"]
    norbitals = params["norbitals"]  # the number of orbitals in the input files
    act_sp = Py2Cpp_int(params["active_space"])
    nsteps = params["nsteps"]
    nfiles = params["nfiles"]
    nstates = len(act_sp)
    istate = params["istate"]
    Nfreqs = params["Nfreqs"]
    ntraj = params["ntraj"]
    set_decoherence = params["set_decoherence"]
    deco_time = params["deco_time"] * units.fs2au
    time_inteval = params["time_inteval"]
    qsh_Ham_prefix = params["qsh_Ham_prefix"]

    rnd = Random()

    ################## Part 1: Collect files and compute statistics ================

    # Electronic Hamiltonian
    H_vib_re = []  # list of MATRIX
    H_vib_im = []  # list of MATRIX
    H_vib = [] # list of CMATRIX

    H_vib_re_ave, H_vib_re_std = None, None
    H_vib_im_ave, H_vib_im_std = None, None
    dw_Hvib_re, up_Hvib_re = None, None
    dw_Hvib_im, up_Hvib_im = None, None

    freqs_re, freqs_im = None, None

    for i in xrange(0, nfiles): # how many files we have

        ##############################################################################
        # Read in the "elementary" overlaps and energies - in the basis of KS orbitals
        ##############################################################################       

        filename_re = params["Hvib_re_prefix"]+str(i)+params["Hvib_re_suffix"]
        filename_im = params["Hvib_im_prefix"]+str(i)+params["Hvib_im_suffix"]
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

    freqs = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]
    for i in xrange(nstates):
        for j in xrange(nstates):
                if i == j:
	            freqs[i][j] =  mat_freqs(H_vib_re, i, j, dt, "out/H_vib_re_E_", Nfreqs) 
	        else:
	            freqs[i][j] = mat_freqs(H_vib_im, i, j, dt, "out/H_vib_im_D_", Nfreqs) 

    
    decoh_times, decoh_rates = decoh_method(set_decoherence, nstates, Hvib, deco_time)    


    if len(freqs[0][0]) < Nfreqs:
        Nfreqs = len(freqs[0][0])
        print "The input Nfreqs is larger than the maximal number of the peaks, changing it to ", Nfreqs

    # Ok, now we have the function - sum of sines, so let's compute the standard deviation
    # This is a silly method - just do it numerically
    dev = [ [ 0.0 for i in xrange(nstates)] for j in xrange(nstates)]
    for r in xrange(1000000):
        fu = [ [ 0.0 for i in xrange(nstates)] for j in xrange(nstates)]
        for i in xrange(nstates):
            for j in xrange(nstates):
                for k in xrange(Nfreqs):
		    fu[i][j] = fu[i][j] + freqs[i][j][k][2] * math.sin(freqs[i][j][k][0]*r*dt/units.au2wavn)
 	
        for i in xrange(nstates):
	     for j in xrange(nstates):
                  dev[i][j] = dev[i][j] + fu[i][j]**2
     
    for i in xrange(nstates):
	 for j in xrange(nstates):
               dev[i][j] = math.sqrt( dev[i][j] / 1000000.0 )
			 
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
    
    # for initization
    Hvib = compute_Hvib(Nfreqs, freqs, t, nstates,
                        H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                        H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                        dev)

    for i in xrange(1, nsteps): # nsteps     

        for traj in xrange(ntraj):
            propagate_electronic(0.5*dt, Cadi[traj], Hvib)
            Cadi[traj] = msdm(Cadi[traj], 0.5*dt, state[traj], decoh_rates)

        # compute QSH Hvib at time i
        Hvib = compute_Hvib(Nfreqs, freqs, i*dt, nstates,
                            H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                            H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                            dev)

        for traj in xrange(ntraj):
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

        # print out the QSH population and vibronc Ham
        out = open("populations.txt", "a")
	out.write("%d " % (i*dt/41.0) )
	for j in xrange(nstates):
	    out.write("   %8.5f  " % pop[j] )
        out.write("\n")
        out.close()

        if i%time_inteval == 0:
	    Hvib.real().show_matrix(qsh_Ham_prefix+str(i)+"_re")
	    Hvib.imag().show_matrix(qsh_Ham_prefix+str(i)+"_im")
      

