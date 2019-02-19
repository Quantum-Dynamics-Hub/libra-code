#*********************************************************************************                     
#* Copyright (C) 2017 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file datautils.py 
# This module implements various utility functions to analyze the arrays of data of different kind
#
# The module contain the following functions:
# find_maxima(s, verbose=0, filename="run.log")
# scalar_stat(X)
# matrix_stat(X)
# matrix_freqs(X, a, b, dt, filename, Nfreqs)

import os
import sys
import math
import copy
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import acf_vector


   


def matrix_freqs(X, a, b, dt, prefix, Nfreqs, verbose = [1,1,1], dw = 1.0, wspan = 3000.0, logfile=None):
# This function analyzes the timeseries of the matrices X
# evolving with the time step dt. The focus will be specifically on the 
# matrix elements with given indices. The function will determine the
# spectrum of the this matrix element. Up to Nfreqs will be determined

# \param[in] X  (list of MATRIX objects) The input data
# \param[in] a, b (int) indices that define the matrix elements to analyze
# \param[in] dt (float) The time step [a.u.] with which the data points in the series 
#            are sampled
# \param[in] prefix (string) A prefix of the filename to which the data will be printed out
# \param[in] Nfreqs (int) The maximal number of frequencies we want to extract
# \param[in] verbose (list of ints) The controls that define which files to print out: 0 - no, 1 - yes
#         verbose[0]:  ACF - format:  <time( in fs)>   <normalized ACF>  <unnormalized ACF>
#         verbose[1]:  FT spectrum - format:  <frequency (in cm^-1)>  <FT amplitude>
#         verbose[2]:  list of ordered frequencies - format: <index> <frequency (in cm^-1)> <FT amplitude> <normalized FT amplitude>
#
# \param[in] dw (float) The spacing between the adjacent points in a computed spectrum [in cm^-1]
# \param[in] wspan (float) The span of frequencies for which the spectrum is to be computed  [in cm^-1]
# \param[in] logfile (string) The name of the file to which we can log some info
#
#
#  Returns: list of ordered FT frequencies and their properties, each entry is a list:
#  [ <index> <frequency (in cm^-1)> <FT amplitude> <normalized FT amplitude> ]

    au2fs = 0.02419 # 40 a.u. is 1 fs 
    inv_cm2ev = (1.0/8065.54468111324)
    au2wavn = 27.211 * 8065.54468111324
    ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
    inv_cm2Ha = inv_cm2ev * ev2Ha
    au2wavn = 27.211 * 8065.54468111324


    dw = dw * inv_cm2Ha            # convert to Ha (atomic units)
    wspan = wspan * inv_cm2Ha      # convert to Ha (atomic units)

    N = len(X)
    sz = X[0].num_of_rows
    freqs = []
            
    # Collect info in a different format
    data_ab = []
    for n in xrange(N):
        data_ab.append(VECTOR(X[n].get(a,b), 0.0, 0.0))
    
    # Now compute ACFs of X matrix elements and print out the corresponding data
    T,  norm_acf,  raw_acf  = acf_vector.acf( acf_vector.center_data(data_ab)  , dt )  # dt is in a.u.


    # Printing the ACF    
    if verbose[0]:
        f = open(prefix+"_acf_"+str(a)+"_"+str(b)+".txt","w")   
        tsz = len(T)
        for it in xrange(tsz):
          f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it]*au2fs , norm_acf[it], raw_acf[it]))
        f.close()
    
    # Do the FT
    W,  J  = acf_vector.ft(norm_acf,  wspan, dw, dt)  # dt is in a.u.
    jsz = len(W)


    # Print the FT spectrum
    if verbose[1]:
        f = open(prefix+"_spectrum_"+str(a)+"_"+str(b)+".txt","w")
        for iw in xrange(jsz):
            f.write("%8.5f  %8.5f  \n" % (W[iw]*au2wavn, J[iw] ) )
        f.close()

    
    # Determine all frequencies of the intensity (peaks) and sort them (in accending manner)
    sp = MATRIX(jsz, 1)
    for iw in xrange(jsz):
        sp.set(iw, J[iw]*J[iw])

    out = find_maxima(sp)


    # Determine the number of actual peaks and use it if it is smaller than the number we requested
    if Nfreqs > len(out):
        Nfreqs = len(out)

    szo = len(out) - 1

    # Compute the sum of the absolute values of the Nfreqs amplitues
    norm = 0.0
    for i in xrange(Nfreqs):
        indx = out[szo-i][0]
        norm = norm + abs(J[indx])

    # Normalize the amplitudes and prepare the output result
    for i in xrange(Nfreqs):
        indx = out[szo-i][0]
        freqs.append( [W[indx]*au2wavn, J[indx], J[indx]/norm ] )


    # Print out this info:
    if verbose[2]:
        f = open(prefix+"_maximal_frequencies_"+str(a)+"_"+str(b)+".txt","w")
        f.write("Maximal peaks from the file "+prefix+"_spectrum_"+str(a)+"_"+str(b)+".txt\n")
       
        for i in xrange(Nfreqs):
            indx = out[szo-i][0]
            f.write("index= %3i  frequency= %8.5f  amplitude= %8.5f normalized_amplitude= %8.5f \n" % (i, W[indx]*au2wavn, J[indx], J[indx]/norm) )
        f.close()    



    if logfile != None:
    
        lgfile = open(logfile, "a")
        print "max frequency for ", prefix, " = ", freqs

        for a in freqs:
            print " Timescale is = ", a[0]/au2wavn, " Ha", " omega = E/hbar ", a[0]/au2wavn, " 2 pi*a.u. of time^-1",\
                  " linear frequency = ", (a[0]/au2wavn)/(2.0*math.pi), " a.u.^-1", " Timescale = ", 2.0*math.pi*au2wavn/a[0], " a.u. of time ",\
                  2.0*math.pi*au2wavn*au2fs/a[0], " fs, Amplitude = ", a[1], " Normalized amplitude = ", a[2]

            lgfile.write("Timescale is = %8.5f Ha, omega = E/hbar %8.5f 2 pi*a.u. of time^-1 \
                          linear frequence = %8.5f a.u.^-1  Timescale = %8.5f a.u. of time \
                          %8.5f fs, Amplitude = %8.5f Normalized amplitude = %8.5f \n" % 
                         ( (a[0]/au2wavn), (a[0]/au2wavn), ((a[0]/au2wavn)/(2.0*math.pi)), (2.0*math.pi*au2wavn/a[0]), (2.0*math.pi*au2wavn*au2fs/a[0]), a[1], a[2] ) 
                        )
        lgfile.close()

   
    return freqs





class TestDatautils(unittest.TestCase):
    def test_1(self):
        """Tests find_maxima(s, verbose=0, filename="run.log") """
        s = MATRIX(6,1)
        slst = [4.0, 20.0, -50.0, 42.0, -2.0]
        for i in xrange(len(slst)):
            s.set(i, slst[i])
          
        out = find_maxima(s)

        self.assertEqual(out, [[1, 20.0], [3, 42.0]])

    def test_2(self):
        """Tests scalar_stat(X) """
        X = [1.0, 2.0, 3.0]
        ave, std = scalar_stat(X)

        self.assertEqual(ave, 2.0)
        self.assertEqual(std, math.sqrt(2.0/3.0))

        X = [1.0, 1.0, 1.0]
        ave, std = scalar_stat(X)

        self.assertEqual(ave, 1.0)
        self.assertEqual(std, 0.0)

    def test_3(self):
        """Tests matrix_stat(X) """

        x00 = [1.0,  2.0, 3.0,  4.0]
        x11 = [1.0, -1.0, 1.0, -1.0]
        x01 = [1.0, -1.0, 2.0, -2.0]
        x10 = [0.0, 1.0,  -3.0, 3.0]

        X = []
        for i in xrange(4):
            m = MATRIX(2,2)
            m.set(0,0, x00[i]); m.set(0,1, x01[i]);
            m.set(1,0, x10[i]); m.set(1,1, x11[i]);
            X.append(m)

        ave, std, dw, up = matrix_stat(X)

        self.assertEqual(ave.get(0,0), 2.5)
        self.assertEqual(ave.get(1,1), 0.0)
        self.assertEqual(ave.get(0,1), 0.0) 
        self.assertEqual(ave.get(1,0), 0.25) 

        self.assertEqual(std.get(0,0), math.sqrt(1.5**2 + 0.5**2 + 0.5**2 + 1.5**2)/2.0)
        self.assertEqual(std.get(1,1), 1.0)
        self.assertEqual(std.get(0,1), math.sqrt(2.5) ) 
        self.assertEqual(std.get(1,0), math.sqrt(0.25**2 + 0.75**2 + 3.25**2 + 2.75**2)/2.0) 

        self.assertEqual(dw.get(0,0),  1.0)
        self.assertEqual(dw.get(1,1), -1.0)
        self.assertEqual(dw.get(0,1), -2.0) 
        self.assertEqual(dw.get(1,0), -3.0) 

        self.assertEqual(up.get(0,0),  4.0)
        self.assertEqual(up.get(1,1),  1.0)
        self.assertEqual(up.get(0,1),  2.0) 
        self.assertEqual(up.get(1,0),  3.0) 


    def test_4(self):
        """Tests matrix_freqs(X, a, b, dt, prefix, Nfreqs, verbose = [1,1,1], dw = 1.0, wspan = 3000.0, logfile=None) """
        pass


        

if __name__=='__main__':
    unittest.main()

