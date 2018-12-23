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


def compute(X, a, b, dt, Nfreqs, filename, logname, dw=1.0, wspan=3000.0):
    """
    Computes the frequencies of the matrix elements

    X          [list of MATRIX] - time-series data
    a,b        [int] - indices that define which matrix elements to analyze
    dt         [double] - time step between the data in the time-series [a.u. of time]
    Nfreqs     [int] - the maximal number of frequencies we want to extract
    filename   [string] - the prefix of the filenames generated
    logname    [string] - the name of the log-file
    dw         [double] - the freqency grinding distance [cm^-1, Default: 1.0 cm^-1]
    wspan      [double] - the width of the spectral window to resolve [cm^-1, Default: 3000 cm^-1]

    Return: [list of lists],
            frequencies:
               freqs[fr][0] - frequency of the mode fr [in cm^-1]
               freqs[fr][1] - amplitude of the mode fr in the influence spectrum [arb. units]
               freqs[fr][2] - normalized amplitued of the mode fr [arb.units]

    """

    nsteps = len(X)    # former N
    sz = X[0].num_of_rows
    freqs = []

    # Collect info in a different format
    data_ab = []
    for n in xrange(nsteps):
        data_ab.append(VECTOR(X[n].get(a,b), 0.0, 0.0))
    
    
    # Now compute ACFs of X matrix elements and print out the corresponding data
#    T,  norm_acf,  raw_acf  = acf_vector.acf( acf_vector.center_data(data_ab)  , dt )  # dt is in a.u.
    T,  norm_acf,  raw_acf  = acf_vector.acf( data_ab  , dt )  # dt is in a.u.
    
    dw = dw * units.inv_cm2Ha           # convert to Ha (atomic units)
    wspan = wspan * units.inv_cm2Ha     # convert to Ha (atomic units)

    
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
    out = comn.find_maxima(sp, logname)

    lgfile = open(logname, "a")
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

    
    lgfile = open(logname, "a")
    for a in freqs:
        lgfile.write(" ========= Mode = %5i =========== \n" % (freqs.index(a)) )
        lgfile.write(" Timescale is = %8.5f [Ha] \n" % (a[0]/units.au2wavn) )
        lgfile.write(" omega = E/hbar = %8.5f [2 pi*a.u. of time^-1] \n" % (a[0]/units.au2wavn)   )
        lgfile.write(" linear frequency = %8.5f [a.u.^-1] \n" % (a[0]/(2.0*math.pi*units.au2wavn))   )
        lgfile.write(" Timescale = %8.5f [a.u. of time] \n " % (2.0*math.pi*units.au2wavn/a[0]) )
        lgfile.write(" Timescale = %8.5f [fs] \n" % ( 2.0*math.pi*units.au2wavn*units.au2fs/a[0]) )
        lgfile.write(" Amplitude = %8.5f \n" % (a[1]) )
        lgfile.write(" Normalized amplitude = %8.5f \n" % (a[2]))

    lgfile.close()

   
    return freqs


