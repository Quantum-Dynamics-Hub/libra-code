#*********************************************************************************
#* Copyright (C) 2018-2019 Wei Li and Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""

.. module:: influence_spectrum
   :platform: Unix, Windows
   :synopsis: This module implements functions to compute the autocorrelation functions 
       and their Fourier spectra (influence spectra) of the time-series of matrices (e.g.
       of the "vibronic" Hamiltonian data sampled along the MD trajectories)

.. moduleauthor:: Wei Li and Alexey V. Akimov


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
import libra_py.units as units
import libra_py.acf_vector as acf_vector
import libra_py.acf_matrix as acf_matrix


def compute_mat_elt(X, a, b, params):
    """Computes the frequencies with which a given matrix element evolves in time

    Args:   
        X ( list of MATRIX ): time-series data
        a ( int ): is the row index of the matrix element to analyze
        b ( int ): is the column index of the matrix element to analyze
        params ( dictionalry ): parameters of the simulation. Contain the following keys:
        
            * **params["dt"]** ( double ): time step between the adjacent data in the time-series [a.u. of time, default: 41.0]
            * **params["nfreqs"]** ( int ): the maximal number of frequencies we want to extract [default: 1]
            * **params["dw"]** ( double ): the freqency grid distance [cm^-1, default: 1.0 cm^-1]
            * **params["wspan"]** ( double ): the width of the spectral window to resolve [cm^-1, default: 4000 cm^-1]
            * **params["do_center"]** ( Boolean ): whether to center the data (True) or not (False) by removing the 
                data-set-averaged values from the current values (then one is dealing with the data fluctuations). [default: True]
            * **params["acf_type"]** ( int ): the definition to compute the ACF:
            
                - 0 : the standard chemical convention [default]
                - 1 : the definition more commonly used in statistics

            * **params["do_output"]** ( Boolean ): select whether to produce files (if True) or not to write any 
                files, but just return the needed values (if False) [default: False]
            * **params["filename"]** ( string ): the prefix of the filenames generated. Doesn't matter 
                if `do_output == False` [default: "influence_spectra_"]
            * **params["logname"]** ( string ): the name of the log-file. Doesn't matter 
                if `do_output == False` [default: "out.log"]             

    Returns:
        list of lists: freqs, where

            * freqs[fr][0] - frequency of the mode fr [in 2*pi*a.u.^-1]
            * freqs[fr][1] - amplitude of the mode fr in the influence spectrum [arb. units]
            * freqs[fr][2] - normalized amplitued of the mode fr [arb.units]

    """

    # Set defaults and check critical parameters
    critical_params = [ ] 
    default_params = { "dt":41.0, "nfreqs":1, "dw":1.0, "wspan":4000.0, "do_center":True, "acf_type":0,
                       "do_output":False, "logname":"out.log", "filename":"influence_spectra_" }
    comn.check_input(params, default_params, critical_params)

 
    # Local variables and dimensions
    nsteps = len(X)    
    sz = X[0].num_of_rows

    dt = params["dt"]
    nfreqs = params["nfreqs"]
    dw = params["dw"] * units.inv_cm2Ha           # convert to Ha (atomic units)
    wspan = params["wspan"] * units.inv_cm2Ha     # convert to Ha (atomic units)    
    do_center = params["do_center"]
    acf_type = params["acf_type"]
    do_output = params["do_output"] 
    filename = params["filename"]
    logname = params["logname"]


    # Collect info in a different format
    data_ab = []
    for n in xrange(nsteps):
        xi = MATRIX(1,1)
        xi.set(0,0, X[n].get(a,b))
        data_ab.append(xi)   
    
    
    #======= Now compute ACFs of X matrix elements and print out the corresponding data =======
    T,  norm_acf,  raw_acf = None, None, None
    data = data_ab
    if do_center:
        data = acf_matrix.center_data(data_ab)
        

    T,  norm_acf,  raw_acf  = acf_matrix.acf( data_ab  , dt, acf_type )  # dt is in a.u.
    
    if do_output:
        f = open(filename+"_acf_"+str(a)+"_"+str(b)+".txt","w")   
        tsz = len(T)
        for it in xrange(tsz):
            f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it]*units.au2fs , norm_acf[it], raw_acf[it]))
        f.close()
    
    #======= Do the FT and print out the corresponding data ==============
    W,  J  = acf_vector.ft(norm_acf,  wspan, dw, dt)  # dt is in a.u.
    jsz = len(W)


    sp = MATRIX(jsz, 1)
    for iw in xrange(jsz):
        sp.set(iw, J[iw]*J[iw])
    
    if do_output:
        f = open(filename+"_spectrum_"+str(a)+"_"+str(b)+".txt","w")
        for iw in xrange(jsz):
            f.write("%8.5f  %8.5f  \n" % (W[iw]*units.au2wavn, J[iw] ) )
        f.close()
    
    #===== Determine all frequencies (peaks) and sort them (in accending manner) ====
    out = comn.find_maxima(sp, logname)

    if do_output:
        lgfile = open(logname, "a")
        lgfile.write("Maximal peaks in the file "+filename+"_spectrum_"+str(a)+"_"+str(b)+".txt\n")


    # Reduce the number of frequencies to the maximal number available
    if nfreqs > len(out):
        nfreqs = len(out)

    szo = len(out) - 1

    #==== Compute the intensities and normalized intensities, do the output =========
    freqs = []
    norm = 0.0
    for i in xrange(nfreqs):
        indx = out[szo-i][0]
        norm = norm + abs(J[indx])
       
    for i in xrange(nfreqs):
        indx = out[szo-i][0]
        freqs.append( [W[indx], J[indx], J[indx]/norm ] )

        if do_output:
            lgfile.write("index= %3i  frequency= %8.5f  amplitude= %8.5f normalized_amplitude= %8.5f \n" % (i, W[indx]*units.au2wavn, J[indx], J[indx]/norm) )
    if do_output:
        lgfile.close()    
    
        lgfile = open(logname, "a")
        for a in freqs:
            lgfile.write(" ========= Mode = %5i =========== \n" % (freqs.index(a)) )
            lgfile.write(" Timescale is = %8.5f [Ha] \n" % (a[0]) )
            lgfile.write(" omega = E/hbar = %8.5f [2 pi*a.u. of time^-1] \n" % (a[0])   )
            lgfile.write(" linear frequency = %8.5f [a.u.^-1] \n" % (a[0]/(2.0*math.pi))   )
            lgfile.write(" Timescale = %8.5f [a.u. of time] \n " % (2.0*math.pi/a[0]) )
            lgfile.write(" Timescale = %8.5f [fs] \n" % ( 2.0*math.pi*units.au2fs/a[0]) )
            lgfile.write(" Amplitude = %8.5f \n" % (a[1]) )
            lgfile.write(" Normalized amplitude = %8.5f \n" % (a[2]))
        lgfile.close()

   
    return freqs, T,  norm_acf,  raw_acf,  W,  J


def compute_all(X, params):
    """Computes the frequencies with which all matrix elements evolve in time

    In particular, we are interested only in the frequencies of the real part of 
    diagonal elements and in frequencies of imaginary part of non-diagonal elements.
    This is a typical situation for the "vibronic" Hamiltonian data in the NA-MD

    Args:
        X ( list of CMATRIX ): time-series data of complex matrices, e.g. "vibronic" Hamiltonian
        params ( dictionary ): parameters controlling the execution of :funct:`compute_mat_elt` function
            ..seealso::`compute_mat_elt` for the full description of the required and allowed parameters and the default values

    Returns:
        list[nstates][nstates][nfreqs][3]: freqs, such that
        
            * freqs[a][b][fr][0] ( double ): frequency of the mode fr for the matrix element X_ab [in 2*pi*a.u.^-1] 
            * freqs[a][b][fr][1] ( double ): amplitude of the mode fr in the influence 
                spectrum for the matrix element X_ab  [arb. units]
            * freqs[a][b][fr][2] ( double ):  normalized amplitued of the mode fr for the matrix 
                element X_ab  [arb.units]

    """

    params_re = dict(params)
    params_re.update({"filename":params["filename"]+"_re_"})

    params_im = dict(params)
    params_im.update({"filename":params["filename"]+"_im_"})
    

    # Split the X array into two arrays - real and imaginary parts
    X_re, X_im = [], []
    nsteps = len(X)
    for step in xrange(nsteps):
        X_re.append(X[step].real())
        X_im.append(X[step].imag())

    # Do the calculations for each matrix element
    nstates = X[0].num_of_cols   

    freqs = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]
    T = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]
    norm_acf = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]
    raw_acf = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]
    W = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]
    J = [ [ [] for i in xrange(nstates)] for j in xrange(nstates)]

 
    for i in xrange(nstates):
        for j in xrange(nstates):
            if i == j:
                freqs[i][j], T[i][j],  norm_acf[i][j],  raw_acf[i][j],  W[i][j], J[i][j] = compute_mat_elt(X_re, i, j, params_re)
            else:
                freqs[i][j], T[i][j],  norm_acf[i][j],  raw_acf[i][j],  W[i][j], J[i][j] = compute_mat_elt(X_im, i, j, params_im)

    return freqs, T,  norm_acf,  raw_acf,  W,  J
