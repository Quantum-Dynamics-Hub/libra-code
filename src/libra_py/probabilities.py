#*********************************************************************************
#* Copyright (C) 2018-2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: probabilities
   :platform: Unix, Windows
   :synopsis: This module aims to implement the probabilities to find system in a given state

.. moduleauthor:: Alexey V. Akimov
  
"""

import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from . import units


def Boltz_quant_prob(E, T):
    """

    Computes the quantum Boltzman probability of occupying different 
    energy levels at given temperature T:  P_i = exp(-E_i/kT) / Z
    where Z = sum_i { exp(-E_i/kT)  }

    Args:
        E ( list of doubles ): energy levels [in a.u.]
        T ( double ): temperature [K]

    Returns:
        double: prob: the probability to  find a system in a discrete state i with
            energy E_i at given temperature T, considering a number of selected energy
            levels as given by the list E
    """
 
    b = 1.0/(units.kB * T)

    nstates = len(E)

    Z = 0.0  # partition function
    prob = []

    for n in range(0,nstates):
        prob.append(math.exp(-E[n]*b))
        Z += prob[n]
    
    for n in range(0,nstates):
        prob[n] = prob[n] / Z

    return prob



def Boltz_cl_prob(E, T):
    """

    Computes the normalized classical Boltzmann probability distribution function 

    Args: 
        E ( double ): the minimum energy level [in a.u.]
        T ( double ): temperature [K]

    Returns:
        double: The probability to have kinetic energy greater than a given threshold value at
            given temperature

    See Also:
        This is essentially a Maxwell-Boltzmann distribution in the energy scale
        Used this: http://mathworld.wolfram.com/MaxwellDistribution.html

    """
 
    x = math.sqrt(E/(units.kB * T))
    res = 1.0

    D = (2.0/(math.sqrt(math.pi) * units.kB * T)) * x * math.exp(-x)
    if D>1.0 or D<0.0:
        print("D = ", D)
        sys.exit(0)
    res = D

    return res




def Boltz_cl_prob_up(E, T):
    """

    Computes the classical Boltzmann probability to have kinetic energy larger than a given 
    threshold E  at temperature T. See Eq. 7 of probabilities_theory.docx.
    The present function is related to it.

    Args: 
        E ( double ): the minimum energy level [in a.u.]
        T ( double ): temperature [K]

    Returns:
        double: The probability to have kinetic energy greater than a given threshold value at
            given temperature

    See Also:
        This is essentially a Maxwell-Boltzmann distribution in the energy scale
        Used this: http://mathworld.wolfram.com/MaxwellDistribution.html

    """
 
    x = math.sqrt(E/(units.kB * T))
    res = 1.0

    D = ERF(x) - math.sqrt(4.0/math.pi) * x * math.exp(-x*x)
    if D>1.0 or D<0.0:
        print("D = ", D)
        sys.exit(0)
    res = 1.0 - D


    return res



def HO_prob(E, qn, T):
    """
    
    Probability that the oscillators are in the given vibrational states
    Multi-oscillator generalization of Eq. 10

    Args:
        E ( list of doubles ): all the energy levels present in the system [in a.u.]
        qn ( list of ints ): quantum numbers for each oscillator
        T ( double ): temperature 

    Returns:
        tuple: (res, prob), where:

            * res ( double ): the probability that the system of N oscillators is in a given
                state, defined by quantum numbers of each oscillator
            * prob ( list of N doubles ): probability with which each of N oscillators occupies 
                given vibrational state
    
    """
    n_freqs = len(E)
 
    res = 1.0
    prob = []
    for i in range(0,n_freqs):
        xi = math.exp(-E[i]/(units.kB*T))        
        prob.append(math.pow(xi, qn[i])*(1.0 - xi))
        res = res * prob[i]

    return res, prob


def HO_prob_up(E, qn, T):
    """

    Probability that the oscillators are in vibrational states with quantum numbers
    above or equal to the minimal quantum numbers provided. Multi-oscillator generalization of Eq. 12

    Args:
        E ( list of doubles ): all the energy levels  present in the system [in a.u.]
        qn ( list of ints ): min quantum numbers for each frequency
        T ( double ): temperature [K]

    Returns:
        tuple: (res, prob), where:

            * res ( double ): the probability that the system of N oscillators is in any state
                higher than given quantum numbers of each oscillators
            * prob ( list of N doubles ): probability with which each oscillator can be found in any of 
                the vibrational states higher than given by the quantum numbers in the qn argument

    """
    n_freqs = len(E)
 
    res = 1.0
    prob = []
    for i in range(0,n_freqs):
        xi = math.exp(-E[i]*qn[i]/(units.kB*T))
        prob.append(xi)
        res = res * prob[i]

    return res, prob


def HO_prob_E_up(E, Emin, T):
    """

    We will compute the probability that a system of N oscillators 
    has energy more or equal of Emin. The oscillators can have only 
    quantized energy values

    Args:
        E ( list of doubles ): all the energy levels  present in the system [in a.u.]
        Emin ( double ): the minimum energy
        T ( double ): temperature [K]

    Note:
        This function is not yet implemented

    """
    pass
    
        
