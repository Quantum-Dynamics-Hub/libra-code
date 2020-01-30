#*********************************************************************************                     
#* Copyright (C) 2017-2019 Brendan Smith, Wei Li, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: acf
   :platform: Unix, Windows
   :synopsis: 
       This module implements the functionality to compute Autocorrelation Functions (ACF)
       and do some transformations of them
       The assumption is that data are provided in a matrix form - not vectors, so we can handle the
       data of arbitrary dimensionality


.. moduleauthor:: Brendan Smith, Wei Li, Alexey V. Akimov

"""

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import units
from . import data_stat


def acf_mat(data, dt, opt=0):
    """Compute the autocorrelation function of the given data set

    Args:
        data ( list of MATRIX(ndof, 1) objects ): sequence of real-valued ndof-dimensional vectors
        dt ( double ): time distance between the adjacent data points [units: general]
        opt ( int ): selector of the convention to to compute ACF

            * 0 : the chemist convention,  (1/(N-h)) Sum_{t=1,N-h} (Y[t]*Y[t+h])
            * 1 : the statistician convention, (1/N) Sum_{t=1,N-h} (Y[t]*Y[t+h])

    Returns:
        tuple: (T, nautocorr, autocorr), where:

            T (list of sz doubles ): lag time scale for the ACF [units: same as for dt]
            nautocorr (list of sz doubles ): normalized ACF
            autocorr (list of sz doubles ): un-normalized ACF

    SeeAlso:
        https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm


    """

    sz = len(data)    # For now, we will use the full data set 

                      ###               how many elements we have in the time series
                      ###  old comments we use only a half of the point, because of the 
                      ###               poorer statistics we get otherwise
    autocorr = []
    ndof = data[0].num_of_rows

    for i in range(0,sz):
        total = 0.0
        for j in range(0,sz-i):
            total += (data[j].T()*data[j+i]).get(0)   # scalar product
        if opt==0:
            autocorr.append( total/((sz-i)*ndof) )  # less bias, chemistry adopted
        elif opt==1:
            autocorr.append( total/(sz*ndof) )      # statistically-preferred option

    #normalize the ACF	
    nautocorr = []
    norm = 1.0
    if math.fabs(autocorr[0])>0.0:
        norm = 1.0/autocorr[0]

    T = []
    for it in range(0,sz):
        T.append(it*dt)
        nautocorr.append( norm * autocorr[it] )

    return T, nautocorr, autocorr




def acf_vec(data, dt, opt=0):
    """Compute the autocorrelation function of the given data set

    Args:
        data ( list of VECTOR objects ): sequence of real-valued 3-dimensional vectors
        dt ( double ): time distance between the adjacent data points [units: general]
        opt ( int ): selector of the convention to to compute ACF

            * 0 : the chemist convention,  (1/(N-h)) Sum_{t=1,N-h} (Y[t]*Y[t+h])
            * 1 : the statistician convention, (1/N) Sum_{t=1,N-h} (Y[t]*Y[t+h])

    Returns:
        tuple: (T, nautocorr, autocorr), where:

            T (list of sz doubles ): lag time scale for the ACF [units: same as for dt]
            nautocorr (list of sz doubles ): normalized ACF
            autocorr (list of sz doubles ): un-normalized ACF

    SeeAlso:
        https://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm


    """

    sz = len(data)    # For now, we will use the full data set 

                      ###               how many elements we have in the time series
                      ###  old comments we use only a half of the point, because of the 
                      ###               poorer statistics we get otherwise
    autocorr = []
    ndof = 3.0

    for i in range(0,sz):
        total = 0.0
        for j in range(0,sz-i):
            total += data[j]*data[j+i]  # scalar product
        if opt==0:
            autocorr.append( total/((sz-i)*ndof) )  # less bias, chemistry adopted
        elif opt==1:
            autocorr.append( total/(sz*ndof) )      # statistically-preferred option

    #normalize the ACF	
    nautocorr = []
    norm = 1.0
    if math.fabs(autocorr[0])>0.0:
        norm = 1.0/autocorr[0]

    T = []
    for it in range(0,sz):
        T.append(it*dt)
        nautocorr.append( norm * autocorr[it] )

    return T, nautocorr, autocorr

