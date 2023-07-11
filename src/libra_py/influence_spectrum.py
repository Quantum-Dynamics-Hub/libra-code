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
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import common_utils as comn
import util.libutil as comn
from . import units
from . import data_stat
from . import acf
from . import ft



def recipe1(data, params):
    """A recipe to compute ACF and its FT for data series
           
    Args:
        data ( list of MATRIX(ndof, 1) objects ): sequence of real-valued ndof-dimensional vectors
        params ( Python dictionary ): controlling the parameters
     
            * **params["dt"]** ( double ): time distance between the adjacent data points [units: fs, default: 1.0]
            * **params["wspan"]** ( double ): window of frequencies for the Fourier transform [ units: cm^-1, default: 3000.0 ]
            * **params["dw"]** ( double ): grid points spacing in the frequency domain [ units: cm^-1, default: 1.0 ]
            * **params["do_output"]** ( Boolean ): whether we print out the data the results into files [ default: False ]
            * **params["acf_filename"]** ( string ): the name of the file where to print the ACF [ default: "acf.txt"]
            * **params["spectrum_filename"]** ( string ): the name of the file where to print the spectrum [ default: "spectrum.txt" ]
            * **params["do_center"]** ( Boolean ): a flag controlling whether to center data (=1) or not (=0)
                Centering means we subtract the average value (over all the data points) from all
                the data points - this way, we convert values into their fluctuations [default: True ]

            * **params["acf_type"]** ( int ): selector of the convention to to compute ACF

                * 0 : the chemist convention,  (1/(N-h)) Sum_{t=1,N-h} (Y[t]*Y[t+h]) [ default ]
                * 1 : the statistician convention, (1/N) Sum_{t=1,N-h} (Y[t]*Y[t+h])

            * **params["data_type"]** ( int ): what is the format of the data?
 
                * 0 : list of MATRIX(ndof, 1) [ default ]
                * 1 : list of VECTOR

    Returns:
        tuple: (T, norm_acf, raw_acf, W, J, J2), where:

            * T ( list of double ): time axis [ units: fs ] 
            * norm_acf ( list of double ): normalized ACF
            * raw_acf ( list of double ): un-normalized ACF
            * W ( list of double ): frequencies axis [ units: cm^-1 ]
            * J ( list of double ): amplitudes of FT 
            * J2 ( list of double ): (1/2pi)*|J|^2 

    """

    critical_params = [ ] 
    default_params = { "dt":1.0, "wspan":3000.0, "dw":1.0, "do_output":False, 
                       "acf_filename":"acf.txt", "spectrum_filename":"spectrum.txt",
                       "do_center":True, "acf_type":0, "data_type":0 }
    comn.check_input(params, default_params, critical_params)


    dt = params["dt"] * units.fs2au            # convert to  atomic units of time
    wspan = params["wspan"] * units.inv_cm2Ha  # convert to Ha (atomic units)
    dw = params["dw"] * units.inv_cm2Ha        # convert to Ha (atomic units)
    do_output = params["do_output"]
    acf_filename = params["acf_filename"]
    spectrum_filename = params["spectrum_filename"]
    do_center = params["do_center"]
    acf_type = params["acf_type"]
    data_type = params["data_type"]
    

    #========
    data_new = data
    if do_center:
        if data_type==0:
            data_new = data_stat.mat_center_data(data)
        elif data_type==1:
            data_new = data_stat.vec_center_data(data)

    #=========== ACFs ==============
    T, norm_acf, raw_acf = None, None, None
    if data_type==0:
        T, norm_acf, raw_acf = acf.acf_mat( data_new , dt, acf_type)
    elif data_type==1:
        T, norm_acf, raw_acf = acf.acf_vec( data_new , dt, acf_type)
    else:
        print("Error: data_type = ", data_type, " is not known\n")
        sys.exit(0)


    sz = len(T)
    for it in range(0,sz):
        T[it] = T[it]/units.fs2au  # convert to fs

    if do_output:
        f = open(acf_filename,"w")
        for it in range(0,sz):
            f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it] , norm_acf[it], raw_acf[it]))
        f.close()


    #=========== FT =============
    W, J = ft.ft(norm_acf, wspan, dw, dt)
    sz = len(W)

    J2 = []
    for iw in range(0,sz):
        W[iw] = W[iw]/units.inv_cm2Ha
        J2.append( (1.0/(2.0*math.pi))*J[iw]*J[iw] )

    if do_output:
        f = open(spectrum_filename,"w")
        for iw in range(0,sz):
            f.write("%8.5f  %8.5f  %8.5f\n" % (W[iw], J[iw], J2[iw] ) )
        f.close()

    return T, norm_acf, raw_acf, W, J, J2



def compute_mat_elt(X, a, b, params):
    """Computes the frequencies with which a given matrix element evolves in time

    Args:   
        X ( list of MATRIX ): time-series data
        a ( int ): is the row index of the matrix element to analyze
        b ( int ): is the column index of the matrix element to analyze
        params ( dictionalry ): parameters of the simulation. Contain the following keys:
            * **params["filename"]** ( string ): the prefix of the filenames generated. Doesn't matter 
                if `do_output == False` [default: "influence_spectra_"]
            * **params["logname"]** ( string ): the name of the log-file. Doesn't matter 
                if `do_output == False` [default: "out.log"]             
            * **params["nfreqs"]** ( int ): the maximal number of frequencies we want to extract [default: 1]

            SeeAlso: recipe1(data, params) for the description of other parameters: 
            
                * dr
                * wspan
                * dw
                * do_output
                * acf_filename
                * spectrum_filename
                * do_center
                * acf_type
                * data_type

    Returns:
        tuple: ( T, norm_acf, raw_acf, W, J, J2, freqs ), where

            freqs ( list of lists ), where

                * freqs[fr][0] - frequency of the mode fr [in 2*pi*a.u.^-1]
                * freqs[fr][1] - amplitude of the mode fr in the influence spectrum [arb. units]
                * freqs[fr][2] - normalized amplitued of the mode fr [arb.units]

            SeeAlso: The first 6 outputs are described in the `recipe1`
                * T
                * norm_acf
                * raw_acf
                * W
                * J
                * J2


    """

    # Set defaults and check critical parameters
    critical_params = [ ] 
    default_params = { "nfreqs":1, "logname":"out.log", "filename":"influence_spectra_", "do_output":0,
                       "do_center":True, "acf_type":1, "data_type":0 }
    comn.check_input(params, default_params, critical_params)

    # Local variables and dimensions
    nfreqs = params["nfreqs"]
    do_output = params["do_output"] 
    filename = params["filename"]
    logname = params["logname"]


    #========= Collect info in a different format =======
    nsteps = len(X)    
    sz = X[0].num_of_rows
    data_ab = []
    for n in range(0,nsteps):
        xi = MATRIX(1,1)
        xi.set(0,0, X[n].get(a,b))
        data_ab.append(xi)   

  
    # ===== Compute the ACFs and their FTs
    params1 = dict(params)
    params1["acf_filename"] = filename+"_acf_"+str(a)+"_"+str(b)+".txt"
    params1["spectrum_filename"] = filename+"_acf_"+str(a)+"_"+str(b)+".txt"
    params1["verbose"] = 0


    T, norm_acf, raw_acf, W, J, J2 = recipe1(data_ab, params1)   # T is in fs, W is in cm^-1
 
    
    #===== Determine all frequencies (peaks) and sort them (in accending manner) ====
    out = data_stat.find_maxima(J2, params1)

    if do_output:
        lgfile = open(logname, "a")
        lgfile.write("Maximal peaks in the file "+params1["spectrum_filename"]+"\n")

    # Reduce the number of frequencies to the maximal number available
    if nfreqs > len(out):
        nfreqs = len(out)

    szo = len(out) - 1

    #==== Compute the intensities and normalized intensities, do the output =========
    freqs = []
    norm, norm2 = 0.0, 0.0
    for i in range(0,nfreqs):
        indx = out[szo-i][0]
        norm = norm + abs(J[indx])
        norm2 = norm2 + J2[indx]
       
    for i in range(0,nfreqs):
        indx = out[szo-i][0]
        freqs.append( [W[indx], J[indx], J[indx]/norm, J2[indx]/norm2 ] )

        if do_output:
            lgfile.write("index= %3i  frequency= %8.5f  amplitude= %8.5f \
                          normalized_amplitude= %8.5f normalized_amplitude2= %8.5f \n" 
                         % (i, W[indx], J[indx], J[indx]/norm, J2[indx]/norm2 ) )
    if do_output:
        lgfile.close()    
    
        lgfile = open(logname, "a")
        for a in freqs:
            lgfile.write(" ========= Mode = %5i =========== \n" % (freqs.index(a)) )
            lgfile.write(" omega = E/hbar = %8.5f [cm^-1] \n" % (a[0])   )
            lgfile.write(" Amplitude = %8.5f \n" % (a[1]) )
            lgfile.write(" Normalized amplitude = %8.5f \n" % (a[2]))
        lgfile.close()

   
    return  freqs, T, norm_acf, raw_acf, W, J, J2



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

    critical_params = [ ] 
    default_params = { "filename":"influence_spectra_" }
    comn.check_input(params, default_params, critical_params)


    params_re = dict(params)
    params_re.update({"filename":params["filename"]+"_re_"})

    params_im = dict(params)
    params_im.update({"filename":params["filename"]+"_im_"})
    

    # Split the X array into two arrays - real and imaginary parts
    X_re, X_im = [], []
    nsteps = len(X)
    for step in range(0,nsteps):
        X_re.append(X[step].real())
        X_im.append(X[step].imag())

    # Do the calculations for each matrix element
    nstates = X[0].num_of_cols   

    freqs = [ [ [] for i in range(0,nstates)] for j in range(0,nstates)]
    T = [ [ [] for i in range(0,nstates)] for j in range(0,nstates)]
    norm_acf = [ [ [] for i in range(0,nstates)] for j in range(0,nstates)]
    raw_acf = [ [ [] for i in range(0,nstates)] for j in range(0,nstates)]
    W = [ [ [] for i in range(0,nstates)] for j in range(0,nstates)]
    J = [ [ [] for i in range(0,nstates)] for j in range(0,nstates)]
    J2 = [ [ [] for i in range(0,nstates)] for j in range(0,nstates)]

 
    for i in range(0,nstates):
        for j in range(0,nstates):
            if i == j:
                freqs[i][j], T[i][j],  norm_acf[i][j],  raw_acf[i][j],  W[i][j], J[i][j], J2[i][j] = compute_mat_elt(X_re, i, j, params_re)
            else:
                freqs[i][j], T[i][j],  norm_acf[i][j],  raw_acf[i][j],  W[i][j], J[i][j], J2[i][j] = compute_mat_elt(X_im, i, j, params_im)

    return freqs, T,  norm_acf,  raw_acf,  W,  J, J2







class TestDatautils(unittest.TestCase):
    def test_1(self):
        """Tests find_maxima(s, verbose=0, filename="run.log") """
        s = MATRIX(6,1)
        slst = [4.0, 20.0, -50.0, 42.0, -2.0]
        for i in range(0,len(slst)):
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
        for i in range(0,4):
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




    
if __name__ == '__main__':

    # Test case: 3 frequences
    data = []
    dt = 1.0 * units.fs2au
    dw = 1.0 * units.inv_cm2Ha
    w1 = 500.0 * units.inv_cm2Ha
    w2 = 1400.0 * units.inv_cm2Ha
    w3 = 850.0 * units.inv_cm2Ha
    wspan = 2000.0 * units.inv_cm2Ha

    for it in range(0,1000):
        t = it * dt
        d = MATRIX(3,1)
        d.set(0, 0, math.sin(w1*t) )
        d.set(1, 0, math.cos(w2*t) )
        d.set(2, 0, math.sin(w3*t) )
        data.append( d )
    
    recipe1(data, 1.0, 2000.0, 1.0)


    # Test case: 3 frequences
    data = []
    for it in range(0,1000):
        t = it * dt
        data.append( VECTOR(math.sin(w1*t), math.cos(w2*t), math.sin(w3*t)) )
    
    recipe2(data, 1.0, 2000.0, 1.0)


    unittest.main()

