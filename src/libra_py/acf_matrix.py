#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#* Copyright (C) 2017 Wei Li, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file acf_matrix.py 
# This module implements the functionality to compute Autocorrelation Functions (ACF)
# and do some transformations of them
# The assumption is that data are provided in a matrix form - not vectors, so we can handle the
# data of arbitrary dimensionality
#
# The module contain the following functions:
# average(data)
# center_data(data)
# center_data2(data)
# acf(data)
# recipe1(data, dt, wspan, dw, acf_filename="acf.txt", spectr_filename="spectrum_.txt", do_center=1)
# recipe2(data, dt, wspan, dw)

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def average(data):
    """    
    This function computes the average value of the data
    # data - (list of MATRIX - (ndof x 1) objects)
    """

    ndof = data[0].num_of_rows  
    ave = MATRIX(ndof, 1)

    sz = len(data)
    for i in xrange(sz):
        ave = ave + data[i] 

    ave = ave/float(sz)
    return ave


def center_data(data):
    """
    This function centers data on zero, by subtracting the average 
    value from each element, dof by dof
    data - (list of MATRIX (ndof x 1) objects) - The data
    """

    data_new = []    
    ave = average(data)
    for d in data:
        data_new.append(d - ave)

    return data_new


def acf(data,dt):
    """
    Compute the autocorrelation function of the given data set

    data - (list of MATRIX (ndof x 1) objects) - Data to analyze
    dt - (float) - time distance between the adjacent data points
    """

    sz = len(data)/2  # how many elements we have in the time series
                      # we use only a half of the point, because of the 
                      # poorer statistics we get otherwise
    autocorr = []
    ndof = data[0].num_of_rows

    for i in range(0,sz):
        total = 0.0
        count = 0.0
        for j in range(0,sz-i):
            total += (data[j].T()*data[j+i]).get(0)   # scalar product
            count += 1.0
        autocorr.append( total/(count*ndof) )


    #normalize the ACF	
    nautocorr = []
    norm = 1.0/autocorr[0]
    T = []
    for it in range(0,sz):
        T.append(it*dt)
        nautocorr.append( norm * autocorr[it] )

    return T, nautocorr, autocorr



def acf2(data,dt):
    """
    Compute the autocorrelation function of the given data set
    using the method described at:
    https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm

    Where C_h = (1/N) Sum_{t=1,N-h} (Y[t]-Y_avg)(Y[t+h]-Y_avg)

    data - (list of MATRIX (ndof x 1) objects) - Data to analyze
    dt - (float) - time distance between the adjacent data points
    """

    N   = len(data)
    avg = MATRIX(data[0].num_of_rows,1)
    for i in xrange(N):
        avg += data[i]
    avg /= float(N)

    T, nC, C  = [], [], []
    #out1 = open("acf.txt","w")
    for k in xrange(N):

        # Compute the autocovariance of data 
        acovar = MATRIX(1,1)
        for i in xrange(N - k):
            acovar += ( (data[i]-avg).T() * (data[i+k]-avg) )
        acovar /= float(N)

        # Compute the variance of data 
        var = MATRIX(1,1)
        for i in xrange(N):
            var += (data[i]-avg).T() * (data[i]-avg)
        var /= float(N)

        C.append(acovar.tr()/var.tr())
        nC.append(C[k]/C[0])
        T.append(k*dt)

    return T, nC, C


def ft(acf_data, wspan, dw, dt):  
    """
    We do have a number of FT and FFT functions in the Libra core, but
    this one may be also convenient to have

    acf_data - (list of floats): C(0), C(dt), C(2*dt), etc. where C - is an ACF
    wspan (float) - is the range of the frequencies we want to compute
    dw (float) - is the distance between the nearby points on the frequency scale
    dt (float) - is the time step
    """

    ############### based on the code from Pyxaid ###################
    sz=len(acf_data)    # the # of input points
    npoints = int(wspan/dw)   # the # of output points    

    J = [0.0] * npoints   # FT
    W = [0.0] * npoints   # frequencies

    for iw in xrange(npoints):
        w = iw * dw

        J[iw] = 1.0  # corresponds to it = 0
        for it in range(1,sz):
            t = it * dt
            J[iw] += 2.0*math.cos(w * t)*acf_data[it]

        W[iw] = w
        J[iw] *= dt

    return W, J


def recipe1(data, dt, wspan, dw, acf_filename="acf.txt", spectrum_filename="spectrum.txt", do_center=1):
    """
    data (MATRIX (ndof x 1) ) - data points (each is a multidimensional)
    dt (float) [ fs ] - timestep between adjacent data points
    dspan (float) [ cm^-1 ] - window of frequencies for the Fourier transform
    dw (float) [ cm^-1 ] - grid points spacing in the frequency domain
    acf_filename (string) - the name of the file where to print the ACF
    spectrum_filename (string) - the name of the file where to print the spectrum
    do_center (int) - a flag controlling whether to center data (=1) or not (=0)
    Centering means we subtract the average value (over all the data points) from all
    the data points - this way, we convert values into their fluctuations.
    """


    # Parameters
    inv_cm2ev = (1.0/8065.54468111324)
    ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
    inv_cm2Ha = inv_cm2ev * ev2Ha
    fs2au = (1.0/0.02419)   # 40 a.u. is 1 fs 

        
    wspan = wspan * inv_cm2Ha  # convert to Ha (atomic units)
    dw = dw * inv_cm2Ha        # convert to Ha (atomic units)
    dt = dt * fs2au            # convert to  atomic units of time

    
    # ACFs
    T, norm_acf, row_acf = acf( center_data(data) , dt)
    sz = len(T)

    f = open(acf_filename,"w")
    for it in xrange(sz):
        f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it]/fs2au , norm_acf[it], row_acf[it]))
    f.close()

    # FT
    W, J = ft(norm_acf, wspan, dw, dt)
    sz = len(W)
    f = open(spectrum_filename,"w")
    for iw in xrange(sz):
        f.write("%8.5f  %8.5f  \n" % (W[iw]/inv_cm2Ha, J[iw] ) )
    f.close()




def recipe2(data, dt, wspan, dw, acf_filename="acf.txt", spectrum_filename="spectrum.txt", do_center=1):
    """
    data (MATRIX (ndof x 1) ) - data points (each is a multidimensional)
    dt (float) [ fs ] - timestep between adjacent data points
    dspan (float) [ cm^-1 ] - window of frequencies for the Fourier transform
    dw (float) [ cm^-1 ] - grid points spacing in the frequency domain
    acf_filename (string) - the name of the file where to print the ACF
    spectrum_filename (string) - the name of the file where to print the spectrum
    do_center (int) - a flag controlling whether to center data (=1) or not (=0)
    Centering means we subtract the average value (over all the data points) from all
    the data points - this way, we convert values into their fluctuations.
    """


    # Parameters
    inv_cm2ev = (1.0/8065.54468111324)
    ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
    inv_cm2Ha = inv_cm2ev * ev2Ha
    fs2au = (1.0/0.02419)   # 40 a.u. is 1 fs 


    wspan = wspan * inv_cm2Ha  # convert to Ha (atomic units)
    dw = dw * inv_cm2Ha        # convert to Ha (atomic units)
    dt = dt * fs2au            # convert to  atomic units of time


    # ACFs
    T, norm_acf, row_acf = acf2(data,dt)
    sz = len(norm_acf)

    f = open(acf_filename,"w")
    for it in xrange(sz):
        f.write("%8.5f  %8.5f  \n" % (T[it]/fs2au, norm_acf[it]))
    f.close()

    # FT
    W, J = ft(norm_acf, wspan, dw, dt)
    sz = len(W)
    f = open(spectrum_filename,"w")
    for iw in xrange(sz):
        f.write("%8.5f  %8.5f  \n" % (W[iw]/inv_cm2Ha, J[iw] ) )
    f.close()


    
if __name__ == '__main__':

    # Parameters
    inv_cm2ev = (1.0/8065.54468111324)
    ev2Ha = (1.0/27.211)    # 27.2 ev is 1 Ha 
    inv_cm2Ha = inv_cm2ev * ev2Ha
    fs2au = (1.0/0.02419)   # 40 a.u. is 1 fs 

    # Test case: 3 frequences
    data = []
    dt = 1.0 * fs2au
    dw = 1.0 * inv_cm2Ha
    w1 = 500.0 * inv_cm2Ha
    w2 = 1400.0 * inv_cm2Ha
    w3 = 850.0 * inv_cm2Ha
    wspan = 2000.0 * inv_cm2Ha

    for it in xrange(1000):
        t = it * dt
        d = MATRIX(3,1)
        d.set(0, 0, math.sin(w1*t) )
        d.set(1, 0, math.cos(w2*t) )
        d.set(2, 0, math.sin(w3*t) )
        data.append( d )
    
    recipe1(data, 1.0, 2000.0, 1.0)
    #recipe2(data, 1.0, 2000.0, 1.0)

