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

"""
  This module aims to implement the probabilities to find system in a given state
"""

import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import units


def Boltz_quant_prob(E, T):
    """
    E [list of doubles] - energy levels [in a.u.]
    T [double] - temperature [K]
 
    Computes the quantum Boltzmann probability of occupying different energy levels at given temperature T

    """
 
    b = 1.0/(units.kB * T)

    nstates = len(E)

    Z = 0.0  # partition function
    prob = []

    for n in xrange(nstates):
        prob.append(math.exp(-E[n]*b))
        Z += prob[n]
    
    for n in xrange(nstates):
        prob[n] = prob[n] / Z

    return prob




def Boltz_cl_prob_up(E, T):
    """
    E [doubles] - the minimum energy level [in a.u.]
    T [double] - temperature [K]

    Computes the classical Boltzmann probability to have kinetic energy larger than a given 
    threshold E  at temperature T

    See Eq. 7. The present function is related to it.
    """
 
    x = math.sqrt(E/(units.kB * T))
    res = 1.0

    """
    This is essentially a Maxwell-Boltzmann distribution in the energy scale
    Used this: http://mathworld.wolfram.com/MaxwellDistribution.html
    """

    D = ERF(x) - math.sqrt(4.0/math.pi) * x * math.exp(-x*x)
    if D>1.0 or D<0.0:
        print "D = ", D
        sys.exit(0)
    res = 1.0 - D


    return res



def HO_prob(E, qn, T):
    """
    E [list of doubles] - all the energy levels present in the system [in a.u.]
    qn [list of ints] - quantum numbers for each frequency
    T [double] - temperature 

    Probability that the oscillators are in the given vibrational states

    Multi-oscillator generalization of Eq. 10
    """
    n_freqs = len(E)
 
    res = 1.0
    prob = []
    for i in xrange(n_freqs):
        xi = math.exp(-E[i]/(units.kB*T))        
        prob.append(math.pow(xi, qn[i])*(1.0 - xi))
        res = res * prob[i]

    return res, prob


def HO_prob_up(E, qn, T):
    """
    E [list of doubles] - all the energy levels  present in the system [in a.u.]
    qn [list of ints]- min quantum numbers for each frequency
    T [double] - temperature [K]

    Probability that the oscillators are in vibrational states with quantum numbers
    above or equal to the minimal quantum numbers provided

    Multi-oscillator generalization of Eq. 12
    """
    n_freqs = len(E)
 
    res = 1.0
    prob = []
    for i in xrange(n_freqs):
        xi = math.exp(-E[i]*qn[i]/(units.kB*T))
        prob.append(xi)
        res = res * prob[i]

    return res, prob


def HO_prob_E_up(E, Emin, T):
    """
    E [list of doubles] - all the energy levels  present in the system [in a.u.]
    Emin [double]- the minimum energy
    T [double] - temperature [K]

    We will compute the probability that a system of N oscillators 
    has energy more or equal of Emin. The oscillators can have only 
    quantized energy values
    """
    pass
    
        
