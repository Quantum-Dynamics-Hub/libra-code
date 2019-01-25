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


import common_utils as comn
import step4
import libra_py.units as units
import libra_py.acf_vector as acf_vector
import influence_spectrum


def compute_freqs(nstates, H_vib_re, H_vib_im, dt, Nfreqs, filename, logname, dw, wspan):
    """
    Compute a matrix of frequencies for each matrix element

    nstates   [int] - the number of states (defines the dimensionality of the matrices)
    H_vib_re  [list of MATRIX] - time-series: real parts of the vibronic Hamiltonian
    H_vib_im  [list of MATRIX] - time-series: imaginary parts of the vibronic Hamiltonian
    dt        [double] - time step between the data in the time-series [a.u. of time]
    Nfreqs    [int] - the maximal number of frequencies we want to extract
    filename  [string] - the prefix of the filenames generated
    logname   [string] - the name of the log-file
    dw        [double] - the freqency grinding distance [cm^-1, Default: 1.0 cm^-1]
    wspan     [double] - the width of the spectral window to resolve [cm^-1, Default: 3000 cm^-1]

    """
    
    freqs = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]
    for i in xrange(nstates):
        for j in xrange(nstates):
            if i == j:
                freqs[i][j] = influence_spectrum.compute(H_vib_re, i, j, dt, Nfreqs, filename+"_re_", logname, dw, wspan)
            else:
                freqs[i][j] = influence_spectrum.compute(H_vib_im, i, j, dt, Nfreqs, filename+"_im_", logname, dw, wspan)

    if len(freqs[0][0]) < Nfreqs:
        Nfreqs = len(freqs[0][0])
        print "The input Nfreqs is larger than the maximal number of the peaks, changing it to ", Nfreqs

    # Ok, now we have the function - sum of sines, so let's compute the standard deviation
    # This is a silly method - just do it numerically
    dev = [ [ 0.0 for i in xrange(nstates)] for j in xrange(nstates)]
        
    for i in xrange(nstates):
        for j in xrange(nstates):

            # Adjust the number of frequencies
            if len(freqs[i][j]) < Nfreqs:
                Nfreqs = len(freqs[i][j])

            fu_ave, fu2_ave = 0.0, 0.0
            for r in xrange(1000000):

                fu = 0.0
                for k in xrange(Nfreqs):
                    fu = fu + freqs[i][j][k][2] * math.sin(freqs[i][j][k][0]*r*dt)

                fu_ave = fu_ave + fu
                fu2_ave = fu2_ave + fu*fu

            dev[i][j] = math.sqrt( (fu2_ave - fu_ave**2) / 1000000.0 ) 	
    
    return freqs, dev


def compute_qs_Hvib(Nfreqs, freqs, t, nstates, 
                 H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re, 
                 H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im, 
                 dev):
    """
    Compute the QSH Hamiltonians

    Nfreqs        [int] - the number of frequencies we want to use in the QSH calculations (upper limit, the actual number could be smaller)
    freqs         [list of lists of doubles] - contains the spectral info for various matrix elements of the sampled Hvib
    t             [double] - time at which we want to reconstruct the QSH [a.u.]
    nstates       [int] - the number of states
    H_vib_re_ave  [MATRIX(nstates, nstates)] - average of energies for direct Hamiltonian
    H_vib_im_ave  [MATRIX(nstates, nstates)] - average of couplings for direct Hamiltonian
    H_vib_re_std  [MATRIX(nstates, nstates)] - std of energies for direct Hamiltonian
    H_vib_im_std  [MATRIX(nstates, nstates)] - std of couplings for direct Hamiltonian
    up_Hvib_re    [MATRIX(nstates, nstates)] - maximum value of energies for direct Hamiltonian
    up_Hvib_im    [MATRIX(nstates, nstates)] - maximum value of couplings for direct Hamiltonian
    dw_Hvib_re    [MATRIX(nstates, nstates)] - minimal value of energies for direct Hamiltonian
    dw_Hvib_im    [MATRIX(nstates, nstates)] - minimal value of couplings for direct Hamiltonian
    std - nstates x nstates matrix, std for direct Hamiltonian
    
    Return: CMATRIX(nstates, nstates) - contains QSH vibronic Hamiltonian at given time

    """
    Hvib_stoch_re = MATRIX(nstates,nstates)
    Hvib_stoch_im = MATRIX(nstates,nstates)

    fu = [ [ 0.0 for i in xrange(nstates)] for j in xrange(nstates)]

    for i in xrange(nstates):
        for j in xrange(nstates):
            # Adjust the number of frequencies
            if len(freqs[i][j]) < Nfreqs:
                Nfreqs = len(freqs[i][j])

            for k in xrange(Nfreqs):
                fu[i][j] = fu[i][j] + freqs[i][j][k][2] * math.sin(freqs[i][j][k][0]*t)


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




def run(H_vib, params):
    """
    The procedure to convert the results of QE/model Hvib calculations to longer timescales using the QSH approach

    H_vib                         [list of lists of CMATRIX] - the vibronic Hamiltonian for all data sets and all time-points

    H_vib[idata][istep].get(i,j) - i,j matrix element for the data set ```idata``` and step in that data set ```istep```

    
    === General control parameters ===

    params["dt"]                  [double] - nuclear dynamics integration time step [in a.u. of time, default: 41.0]
    params["nfreqs"]              [int] - maximal number of frequencies to use to reconstruct all the matrix elements
    params["dw"]                  [double] - frequency spacing [cm^-1, Default: 1.0]
    params["wspan"]               [double] - maximal frequency value [cm^-1, Default: 4000.0]
    params["filename"]            [string] - prefix for the filenames where the spectral calculations will be printed out
    params["logname"]             [string] - the name of the file to contain general output of the spectral calculations

    === Required by the input ===
   
    params["nstates"]             [int] - how many lines/columns in the file - the total number of states
    params["active_space"]        [list of ints] - the indices of the states we care about, default: range(nstates)
    params["nfiles"]              [int] - how many input files (direct Hvib) to read, starting from index 0
    params["data_set_paths"]      [string] - where the input files are located
    params["Hvib_re_prefix"]      [string] - prefixes of the files with real part of the vibronic Hamiltonian at time t
    params["Hvib_re_suffix"]      [string] - suffixes of the files with real part of the vibronic Hamiltonian at time t
    params["Hvib_im_prefix"]      [string] - prefixes of the files with imaginary part of the vibronic Hamiltonian at time t
    params["Hvib_im_suffix"]      [string] - suffixes of the files with imaginary part of the vibronic Hamiltonian at time t

    === Define the output ===
    params["do_output"]           [Boolean] - the flag that determines whether to generate the output files. Default: True. 
                                  Set to False for the on-the-fly QSH
    params["nsteps"]              [int] - how many output files (QSH) to write, starting from index 0
    params["output_set_paths"]    [string] - where the resulting Hvib files are to be stored [default: the same as the input paths]
    params["qsh_Hvib_re_prefix"]  [string] - prefixes of the output files with real part of the QSH vibronic Hamiltonian at time t
    params["qsh_Hvib_re_suffix"]  [string] - suffixes of the output files with real part of the QSH vibronic Hamiltonian at time t
    params["qsh_Hvib_im_prefix"]  [string] - prefixes of the output files with imaginary part of the QSH vibronic Hamiltonian at time t
    params["qsh_Hvib_im_suffix"]  [string] - suffixes of the output files with imaginary part of the QSH vibronic Hamiltonian at time t

    Return: list of reconstructed QSH objects for 
     
    """

    critical_params = [ "nstates", "nfiles", "nsteps"] # "data_set_paths", "output_set_paths" ]
    default_params = { "dt":41.0, "nfreqs":1, "dw":1.0, "wspan":4000.0, "logname":"out.log",
                       "filename":"influence_spectra_",
                       "Hvib_re_prefix":"Hvib_", "Hvib_im_prefix":"Hvib_",
                       "Hvib_re_suffix":"_re",   "Hvib_im_suffix":"_im",
                       "qsh_Hvib_re_prefix":"qsh_Hvib_", "qsh_Hvib_im_prefix":"qsh_Hvib_",
                       "qsh_Hvib_re_suffix":"_re",   "qsh_Hvib_im_suffix":"_im",
                       "do_output":True                 }
    comn.check_input(params, default_params, critical_params)

    """ 
    if(len(params["data_set_paths"]) != len(params["output_set_paths"])):
        print "Error: Input and output sets paths should have equal number of entries\n"
        print "len(params[\"data_set_paths\"]) = ", len(params["data_set_paths"])
        print "len(params[\"output_set_paths\"]) = ", len(params["output_set_paths"])
        print "Exiting...\n"
        sys.exit(0)
    """

    dt = params["dt"]
    nfreqs = params["nfreqs"]

    filename = params["filename"] 
    logname = params["logname"]
    dw = params["dw"]
    wspan = params["wspan"]

                             
    nfiles = params["nfiles"]
    nsteps = params["nsteps"]
    nstates = params["nstates"]
    ndata = len(params["data_set_paths"])
    ntraj = params["ntraj"]
    nitimes = len(params["init_times"])
    Ntraj = ndata * nitimes * ntraj

    


    qsh_H_vib = []

    for idata in xrange(ndata):   # over all MD trajectories (data sets)
        
        #======== Read in the vibronic Hamiltonian along the trajectory for each data set ============

        H_vib_re = []  # list of MATRIX
        H_vib_im = []  # list of MATRIX

        for i in xrange(nfiles):
            H_vib_re.append(H_vib[idata][i].real())
            H_vib_im.append(H_vib[idata][i].imag())


        #======== Analyze the Hvib time-seris  ============
        H_vib_re_ave, H_vib_re_std, dw_Hvib_re, up_Hvib_re = comn.mat_stat(H_vib_re)
        H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im = comn.mat_stat(H_vib_im)    

        freqs, dev = compute_freqs(nstates, H_vib_re, H_vib_im, dt, nfreqs, filename, logname, dw, wspan)

        
        #============= Output the resulting QSH Hamiltonians ===========================
        Hvib = []
        for i in xrange(nsteps):
            # compute QSH Hvib at time t_i = i * dt
            qs_Hvib = compute_qs_Hvib(nfreqs, freqs, i*dt, nstates,  H_vib_re_ave, H_vib_re_std, dw_Hvib_re, 
                            up_Hvib_re, H_vib_im_ave, H_vib_im_std, dw_Hvib_im, up_Hvib_im,  dev)

            Hvib.append(CMATRIX(qs_Hvib))

            if params["do_output"]==True:
                #============= Output the resulting QSH Hamiltonians ===========================
                re_filename = prms["output_set_paths"][idata] + prms["qsh_Hvib_re_prefix"] + str(i) + prms["qsh_Hvib_re_suffix"]
                im_filename = prms["output_set_paths"][idata] + prms["qsh_Hvib_im_prefix"] + str(i) + prms["qsh_Hvib_im_suffix"]        
                qs_Hvib.real().show_matrix(re_filename)
                qs_Hvib.imag().show_matrix(im_filename)

        qsh_H_vib.append(Hvib)        
        

    return qsh_H_vib

