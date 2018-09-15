#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
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
import unittest


# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



def hermite(n, x):

    r,s,t = 0.0, 0.0, 0.0
    p,q = 1.0, 0.0

    for m in xrange(n):
        r,s = p,q
        p = 2.0*x*r - 2.0*m*t
        q = 2.0*(m+1)*r
        t = r

    return p


def ket_n(q, n, k, m):
    """
    HO state |n>
    """

    hbar = 1.0  # atomic units
    omega = math.sqrt(k/m)  
    alp = m*omega/hbar

    N_n =  math.pow(alp/math.pi, 0.25) / math.sqrt(math.pow(2.0, n) * FACTORIAL(n))
    ksi = math.sqrt(alp)*q    
    H_n = hermite(n, ksi)
 
    res = N_n * H_n * math.exp(-0.5*ksi*ksi)

    return res
        




def HO_sup(q, params):
    """
    The probability density function: superposition of HO eigenstates

    """

    k = params["k"]
    m = params["m"]
    states = params["states"]
    coeffs = params["coeffs"]

    x = q.get(0)

    sz = len(states)
    p = 0.0
    for n in xrange(sz):
        p = p + coeffs[n] * ket_n(x, states[n], k, m)

    p = p * p 

    return p


def HO_sup_t(q, params, t):
    """
    The probability density function: superposition of HO eigenstates

    now, time-dependent

    |Psi> = \sum_n {  c_n * |n> * exp(-i*t*E_n/hbar) }

    """

    k = params["k"]
    m = params["m"]
    states = params["states"]
    coeffs = params["coeffs"]

    hbar = 1.0  # atomic units
    omega = math.sqrt(k/m)  


    x = q.get(0)

    sz = len(states)
    p = 0.0+0.0j
    for n in xrange(sz):
        E_n = hbar*omega*(n+0.5)
        p = p + coeffs[n] * ket_n(x, states[n], k, m) * cmath.exp(-1.0j*t*E_n/hbar)

    p = (p.conjugate() * p ).real

    return p



def bin(sample, min_, max_, dx, i, j, filename):

    # Prepare the grids
    x_points, y_points = [], []
    max_pts = int((max_ - min_)/dx) + 1

    for n in xrange(max_pts):
        x_points.append(min_ + n * dx)
        y_points.append(0.0)

    # Compute the frequencies
    sz = len(sample)
    for n in xrange(sz):
        x = sample[n].get(i,j) 
        indx = int((x - min_)/dx)

        y_points[indx] = y_points[indx] + 1.0/float(sz)
              

    f = open(filename, "w")
    for n in xrange(max_pts):
        f.write("%8.5f  %8.5f \n" % (x_points[n], y_points[n]))
    f.close()


def test_1():

    rnd = Random()

    q = MATRIX(1,1);  q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[0], "coeffs":[1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-1.txt")    

    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[1], "coeffs":[1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-2.txt")    

    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[2], "coeffs":[1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-3.txt")    


    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[5], "coeffs":[1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-4.txt")    

    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[10], "coeffs":[1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-5.txt")    

    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[25], "coeffs":[1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-6.txt")    


    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[0, 1], "coeffs":[1.0, 1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-7.txt")    


    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[0, 1], "coeffs":[1.0, -1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-8.txt")    

    q.set(0, 0.5)    
    params = {"k":1.0, "m":2000.0, "states":[0, 1, 2], "coeffs":[1.0, 1.0, 1.0]}    
    sampling = metropolis_gau(rnd, HO_sup, q, params, 1000000, 50000, 0.05) 
    print len(sampling)
    bin(sampling, -1.5, 2.0, 0.01, 0, 0, "_distrib-9.txt")    


           
test_1()






