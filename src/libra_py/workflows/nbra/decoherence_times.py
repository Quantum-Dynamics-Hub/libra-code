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

"""

  This module takes care of the decoherence times calculations

"""

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import units
import common_utils as comn


__all__ = ["energy_gaps",
           "energy_gaps_ave", 
           "decoherence_times",
           "decoherence_times_ave"
          ]


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


def energy_gaps_ave(Hvib, itimes, nsteps):
    """
    Pre-compute the energy gaps along the trajectory 

    Hvib          [list of lists of CMATRIX] - Vibronic Hamiltonians along the trajectory for different data sets (adiabatic MDs)
    itimes        [list if ints]             - initial times for averaging
    nsteps        [int]                      - the length of the sub-data

    """
    
    ndata = len(Hvib)
    nitimes = len(itimes)
    nstates = Hvib[0].num_of_cols
    
    dE = []
    for step in xrange(0, nsteps):
        dEij = MATRIX(nstates, nstates)

        for i in xrange(nstates):
            for j in xrange(i+1, nstates):

                de = 0.0
                for idata in xrange(ndata):
                    for it_indx in xrange(nitimes): 
                        it = itimes[it_indx]

                        deij = deij + math.fabs(Hvib[idata][it+step].get(i,i).real - Hvib[idata][it+step].get(j,j).real)
                de = de/float(nitimes*ndata)

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
    dE_ave, dE_std, dE_dw_bound, dE_up_bound = comn.mat_stat(dE)

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





def decoherence_times_ave(Hvib, itimes, nsteps, verbosity=0):
    """
    Hvib          [list of lists of CMATRIX] - Vibronic Hamiltonians along the trajectory for different data sets (adiabatic MDs)
    itimes        [list if ints]             - initial times for averaging
    nsteps        [int]                      - the length of the sub-data

    Compute the decoherence times:
    Ref: Akimov, A. V; Prezhdo O. V. J. Phys. Chem. Lett. 2013, 4, 3857  

    """

    # Compute energy gaps
    dE = energy_gaps_ave(Hvib, itimes, nsteps)
    dE_ave, dE_std, dE_dw_bound, dE_up_bound = comn.mat_stat(dE)

    nstates = Hvib[0][0].num_of_cols
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


