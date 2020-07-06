#*********************************************************************************
#* Copyright (C) 2020 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: gaussian_kernel_algorithm
   :platform: Unix, Windows
   :synopsis: This module implements functions for the Gaussian kernel algorithm 
       discussed in the paper: arXiv:1712.01918 [physics, physics:quant-ph]

.. moduleauthor:: 
       Alexey V. Akimov 
  
"""


import os
import sys
import math
import cmath
import re

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn


def compute_apriory_prob_densities_1D(q):
    """
    This is a function to compute the probability density of the set of 1D points
    via Voronoi tesselattion. 
    
    Args:
    
        q (npts x float): coordinates of each 1D trajectory
        
    Returns:    
        list: prob_dens        
            * prob_dens (npts x float): the probability densities corresponding to each input point,
                in the same order as the original data
             
    References:
        arXiv:1712.01918 [physics, physics:quant-ph], alternative to Eq. 31
        
    """
    
    npts = len(q)    
                        
    proj = []        
    for pt in range(npts):
        proj.append( [pt, q[pt] ] )            
        sorted_proj = merge_sort(proj)
    
        
    prob_dens = []             # along the original input

    for ipt in range(npts):
        prob_dens.append(0.0)

            
    for ipt in range(npts):                    
        pt = sorted_proj[ipt][0]       # the original index of the point
        
        if ipt==0:            
            x = sorted_proj[ipt][1]        # the coordinate of this point       
            x_next = sorted_proj[ipt+1][1] # the coordinate of the next point                    
            dx = (x_next - x)                         
        elif ipt==npts-1:     
            x_prev = sorted_proj[ipt-1][1] # the coordinate of the next point                    
            x = sorted_proj[ipt][1]        # the coordinate of this point       
            dx = (x - x_prev)             
        else:            
            x_prev = sorted_proj[ipt-1][1] # the coordinate of the previous point
            x_next = sorted_proj[ipt+1][1] # the coordinate of the next point        
            dx = 0.5 * (x_next - x_prev)


        prob_dens[pt] = 1.0/((npts+1)*dx) 
        
                    
    return prob_dens
    


def gaussian_density_kernel_1D(x):
    """
    This function is a simple Gaussian kernel for probability density. 
    
    Args:
        x (double): input value, 1D
        
    Returns:
        double: the value of the normalized 1D Gaussian function at the input argument value
        
    References:
        arXiv:1712.01918 [physics, physics:quant-ph], after Eq. 27     
    
    """    
    res = (1.0/math.sqrt(2.0*math.pi) ) * math.exp(-0.5*x*x)
    
    return res


def gaussian_density_estimator_1D(x, q, h):
    """
    This function computes the estimate of the probability density
    at a given point `X` as created by the superposition of multiple 
    1D Gaussian density kernels, centeres at the points `Q` and having the
    width parameters `h`.
    
    Args:
        x (double): coordinate of the point of interest
        q (npts x double): the centers of the 1D Gaussian kernels ("worlds")
        h (npts x double): 1D Gaussian width parameters for each center
        
    Returns: 
        double: the value of the probability density at a given point
        
    References:
        arXiv:1712.01918 [physics, physics:quant-ph], Eq. 27
    """
    
    npts = len(q) # number of points
    
    res = 0.0
    for i in range(npts):
        res += gaussian_density_kernel_1D( (x-q[i])/h[i] ) / h[i]
        
    return res / npts

        
def gaussian_kernel_algorithm_iteration_1D(q, h, target_probs):
    """
    This function refines the input Gaussian kernel width parameters, `h`,
    to better match the probability density values computed by the estimator
    at given set of points, `q`, to the target probability densities `target_probs`.
    In this implementation, the probability density estimator is also defined in
    terms of all points in `q`.
    
    Args:
        q (npts x double): the centers of the 1D Gaussian kernels ("worlds")
        h (npts x double): 1D Gaussian width parameters for each center
        target_probs (npts x double): values of the target (apriory probability) at `q` points
    
    Returns:
        None: but the algorithm updates the input `h` parameters
    
    References:
        arXiv:1712.01918 [physics, physics:quant-ph], Eq. 32
    """
    
    npts = len(q)
    
    h0 = list(h)
    
    for i in range(npts):
        h[i] = h0[i] * gaussian_density_estimator_1D(q[i], q, h0) / target_probs[i]
        


def compute_widths_1D(q, niter=5, guess_h_val = 1.0):
    """
    This function computes the best set of width parameters for the Gaussian
    kernel functions that approximate the probability density of a given distribution of the points
    
    Args:
        q (npts x double): the centers of the 1D points, whose probability density function we want compute
        niter (int) : the number of refinement iterations
        guess_h_val (double): the guess value of the width parameters
    
    Returns:
        list: h
            * h (npts x doubles): the width parameters for each of the 1D Gaussian kernel functions, ordered
                according to the order of the original input centers 

        
    References:
        arXiv:1712.01918 [physics, physics:quant-ph], Eq. 32        
    """
    npts = len(q)
               
    target_dens = compute_apriory_prob_densities_1D(q)
    
    h = []
    for i in range(npts):
        h.append(guess_h_val)
        
    # Iterations
    for i in range(niter):
        gaussian_kernel_algorithm_iteration_1D(q, h, target_dens)
        
    return h
        

    