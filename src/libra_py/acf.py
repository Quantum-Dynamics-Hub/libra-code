#*********************************************************************************                     
#* Copyright (C) 2017 Wei Li, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file acf.py 
# This module implements the functionality to compute Autocorrelation Functions (ACF)
# and do some transformations of them
#
# The module contain the following functions:
# average(data)
# center_data(data)
# acf(data)

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def average(data):
## This function computes the average value of the data
# data - is a list of VECTOR objects
    ave = VECTOR()
    sz = len(data)

    for i in xrange(sz):
        ave = ave + data[i] 

    ave = ave/float(sz)
    return ave


def center_data(data):
## this function centers data on zero, by subtracting the average 
# value from each element
# data - is a list of VECTOR objects

    data_new = []    
    ave = average(data)
    for d in data:
        data_new.append(d - ave)

    return data_new


def acf(data,dt):
## data - is a list of VECTOR objects

    sz = len(data)/2  # how many elements we have in the time series
                      # we use only a half of the point, because of the 
                      # poorer statistics we get otherwise
    autocorr = []

    for i in range(0,sz):
        total = 0.0
        count = 0.0
        for j in range(0,sz-i):
            total += data[j]*data[j+i] 
            count += 1.0
        autocorr.append( total/count )


    #normalize the ACF	
    nautocorr = []
    norm = 1.0/autocorr[0]
    T = []
    for it in range(0,sz):
        T.append(it*dt)
        nautocorr.append( norm * autocorr[it] )

    return T, nautocorr, autocorr


def ft(acf_data, wspan, dw, dt):  
# we do have a number of FT and FFT functions in the Libra core, but
# this one may be also convenient to have
# acf_data - list of floats: C(0), C(dt), C(2*dt), etc. where C - is an ACF
# span - is the range of the frequencies we want to compute
# dw - is the distance between the nearby points on the frequency scale
# dt - is the time step


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


def recipe1(data, dt, wspan, dw):
# dt in fs
# dspan in cm^-1
# dw in cm^-1


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

    f = open("acf.txt","w")
    for it in xrange(sz):
        f.write("%8.5f  %8.5f  %8.5f  \n" % (T[it]/fs2au , norm_acf[it], row_acf[it]))
    f.close()

    # FT
    W, J = ft(norm_acf, wspan, dw, dt)
    sz = len(W)
    f = open("spectrum_.txt","w")
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
        data.append( VECTOR(math.sin(w1*t), math.cos(w2*t), math.sin(w3*t)) )
    
    recipe1(data, 1.0, 2000.0, 1.0)
