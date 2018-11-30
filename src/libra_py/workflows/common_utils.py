#*********************************************************************************
#* Copyright (C) 2017-2018 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  
#
import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import units

    
def get_matrix(nrows, ncols, filename_re, filename_im, act_sp):
    """
    This file reads the real and imaginary components of a matrix of given original size,  
    takes its sub-matrix (as defined by the act_sp function) and returns the resulting 
    complex matrix

    nrows (int)   - the number of rows in the original matrix (read from the files)
    ncols (int)   - the number of columns in the original matrix (read from the files)
    filename_re (string) - the name of the file containing the real part of the matrix 
    filename_im (string) - the name of the file containing the imaginary part of the matrix 
    act_sp (list of ints) - the indices of the columns and rows to be taken to construct the resulting matrices
    
    """

    X_re = MATRIX(nrows, ncols); X_re.Load_Matrix_From_File(filename_re)
    X_im = MATRIX(nrows, ncols); X_im.Load_Matrix_From_File(filename_im)

    nstates = len(act_sp)
    x_re = MATRIX(nstates, nstates);
    x_im = MATRIX(nstates, nstates);

    pop_submatrix(X_re, x_re, act_sp, act_sp)
    pop_submatrix(X_im, x_im, act_sp, act_sp)

    return CMATRIX(x_re, x_im)





def find_maxima(s):
    """
    s [list of double] - data

    This function finds all the maxima of the data set and sorts them according to the data
    The maxima are defined as s[i-1] < s[i] > s[i+1]
   
    Returns the list of the indices of the maximal values
    """

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
    """
    Computes the average and std 
    X [list of double] - data
    """
 
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
    """
    Computes the average and std     
    X [list of MATRIX] - data
    """

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


def cmat_stat(X):
    """
    Computes the average and std     
    X [list of CMATRIX] - data
    """

    N = len(X)
    res = CMATRIX(X[0]); res *= 0.0

    #===== Average ====
    for i in xrange(N):
        res = res + X[i]
    res = res / float(N)


    #===== Std ========
    res2 = CMATRIX(X[0]); res2 *= 0.0

    for a in xrange(res2.num_of_rows):
        for b in xrange(res2.num_of_cols):
        
            tmp = 0.0+0.0j
            for i in xrange(N):
                dx = X[i].get(a,b) - res.get(a,b)
                tmp = tmp + (dx.conjugate() * dx)
            tmp = math.sqrt( tmp / float(N) )

            res2.set(a,b, tmp)

    return res, res2



def energy_gaps(Hvib):
    """
    Pre-compute the energy gaps along the trajectory 

    Hvib [list of CMATRIX] - Vibronic Hamiltonians along the trajectory
    """

    nsteps = len(Hvib)
    nstates = Hvib[0].num_of_cols
    
    dE = []
    for step in xrange(0, nsteps):
        dEij = MATRIX(nstates, nstates)

        for i in xrange(nstates):
            for j in xrange(i+1, nstates):

                deij = math.fabs(Hvib[step].get(i,i).real - Hvib[step].get(j,j).real)
                dEij.set(i,j, deij)
                dEij.set(j,i, deij)

        dE.append(dEij)

    return dE




def decoherence_times(Hvib, verbosity=0):
    """
    Hvib [list of CMATRIX] - timeseries of the vibronic Hamiltonian

    Compute the decoherence times:
    Ref: Akimov, A. V; Prezhdo O. V. J. Phys. Chem. Lett. 2013, 4, 3857  

    """

    # Compute energy gaps
    dE = energy_gaps(Hvib)
    dE_ave, dE_std, dE_dw_bound, dE_up_bound = mat_stat(dE)

    nstates = Hvib[0].num_of_cols
    decoh_times = MATRIX(nstates, nstates)
    decoh_rates = MATRIX(nstates, nstates)

    for a in xrange(nstates):
        for b in xrange(nstates):
            if a==b:
                decoh_times.set(a,a, 1000000.0)
                decoh_rates.set(a,a, 0.0)
            else:
                de = dE_std.get(a,b)
                if de>0.0:
                      tau = math.sqrt(12.0/5.0) / de
                      decoh_times.set(a,b, tau)
                      decoh_rates.set(a,b, 1.0/tau)

    if verbosity>0:
        print "Decoherence times matrix (a.u. of time):"
        decoh_times.show_matrix()

        print "Decoherence times matrix (fs):"
        tmp = decoh_times * units.au2fs
        tmp.show_matrix()

        print "Decoherence rates matrix (a.u.^-1):"
        decoh_rates.show_matrix()


    return decoh_times, decoh_rates

