#*********************************************************************************
#* Copyright (C) 2017-2019 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
.. module:: decoherence_times
   :platform: Unix, Windows
   :synopsis: 
       This module implements functions to compute decoherence times and relevant quantities

.. moduleauthor:: Alexey V. Akimov

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
from libra_py import data_stat

__all__ = ["decoherence_times2rates",
           "energy_gaps",
           "energy_gaps_ave", 
           "decoherence_times",
           "decoherence_times_ave"
          ]



def decoherence_times2rates(tau):
    """

    An auxiliary function to convert the decoherence times matrix into the 
    decoherence rates matrix
     
    Args:
        tau ( MATRIX(N,N) ): the matrix of decoherence times, diagonal elements are
                set to very large number which corresponds to no decoherence of a state with itself [ units: a.u. ]

    Returns:
        MATRIX(N,N): decoh_rates:  the matrix of decoherence rates, diagonal elements are
            set to zero which corresponds to no decoherence of a state with itself, the off-diagonal
             elements are equal to inverse of the off-diagonal matrix elements of ```tau``` [ units: a.u.^-1 ]
    """

    nstates = tau.num_of_cols
    decoh_rates = MATRIX(nstates, nstates)

    for i in range(0,nstates):
        for j in range(0,nstates):
            if i==j:
                decoh_rates.set(i,i, 0.0)
            else:
                tau_ij = tau.get(i,j)
                if tau_ij > 0.0:                      
                    decoh_rates.set(i,j, 1.0/tau_ij)

    return decoh_rates



def energy_gaps(Hvib):
    """Pre-compute the energy gaps along the trajectory 

    Args:     
        Hvib ( list of MATRIX objects ): Vibronic Hamiltonians along the trajectory

    Returns:
        ( list of MATRIX(nstates, nstates) ): dE, where:
            dE[t].get(i,j) is the energy gap between states i and j at time t:
            E_i(t) - E_j(t) 

    """

    nsteps = len(Hvib)
    nstates = Hvib[0].num_of_cols
    
    dE = []
    for step in range(0, nsteps):
        dEij = MATRIX(nstates, nstates)

        for i in range(0,nstates):
            for j in range(i+1, nstates):

                deij = math.fabs(Hvib[step].get(i,i).real - Hvib[step].get(j,j).real)
                dEij.set(i,j, deij)
                dEij.set(j,i, deij)

        dE.append(dEij)

    return dE


def energy_gaps_ave(Hvib, itimes, nsteps):
    """Pre-compute the energy gaps along the trajectory 

    Args:
        Hvib ( list of lists of CMATRIX objects ):
            Vibronic Hamiltonians along the trajectory for different data sets (adiabatic MDs),
            and potential different sections of the time-range
            where Hvib[idata][istep] is a CMATRIX object that represents a vibronic Hamiltonian
            from the data set ```idata``` at the time step ```istep```. 

        itimes ( list if ints ): initial times for averaging. The number ```itimes[idata]``` tells 
            which datapoint (timestep) of the time-series Hvib[idata] consider the beginning of 
            the range that will be used to compute gap fluctuations

        nsteps ( int ): the length of the time-steps to consider in the calculations for each data set

    Returns:
        ( list of MATRIX(nstates, nstates) ): dE, where:
            dE[t].get(i,j) is the absolute value of energy gap between states i and j at time t, but also averaged 
            over several data sets:      < |E_i(t) - E_j(t)| >

    """
    
    ndata = len(Hvib)
    nitimes = len(itimes)
    nstates = Hvib[0][0].num_of_cols
    
    dE = []
    for step in range(0, nsteps):
        dEij = MATRIX(nstates, nstates)

        for i in range(0,nstates):
            for j in range(i+1, nstates):

                deij = 0.0
                for idata in range(0, ndata):
                    for it_indx in range(0, nitimes): 
                        it = itimes[it_indx]

                        deij = deij + math.fabs(Hvib[idata][it+step].get(i,i).real - Hvib[idata][it+step].get(j,j).real)
                deij = deij/float(nitimes*ndata)

                dEij.set(i,j, deij)
                dEij.set(j,i, deij)

        dE.append(dEij)

    return dE




def decoherence_times(Hvib, verbosity=0):
    """Compute the matrix of decoherence times from the time-series data

    Ref: Akimov, A. V; Prezhdo O. V. J. Phys. Chem. Lett. 2013, 4, 3857  


    Args:
        Hvib ( list of CMATRIX objects ): timeseries of the vibronic Hamiltonian
        verbosity ( int ): the flag controlling the amount of extra output
            Value of 0 [ default ] prints no additional output

    Returns:
        tuple:  ( decoh_times, decoh_rates ), where

            * decoh_times ( MATRIX(N,N) ): the matrix of decoherence times, diagonal elements are
                set to very large number which corresponds to no decoherence of a state with itself [ units: a.u. ]

            * decoh_rates ( MATRIX(N,N) ): the matrix of decoherence rates, diagonal elements are
                set to zero which corresponds to no decoherence of a state with itself, the off-diagonal
                elements are equal to inverse of the off-diagonal matrix elements of ```decoh_times``` [ units: a.u.^-1 ]

    """

    # Compute energy gaps
    dE = energy_gaps(Hvib)
    dE_ave, dE_std, dE_dw_bound, dE_up_bound = data_stat.mat_stat(dE)

    nstates = Hvib[0].num_of_cols
    decoh_times = MATRIX(nstates, nstates)

    for a in range(0,nstates):
        for b in range(0,nstates):
            if a==b:
                decoh_times.set(a,a, 1.0e+10)
            else:
                de = dE_std.get(a,b)
                if de>0.0:
                      tau = math.sqrt(12.0/5.0) / de
                      decoh_times.set(a,b, tau)

    decoh_rates = decoherence_times2rates(decoh_times)

    if verbosity>0:
        print("Decoherence times matrix (a.u. of time):")
        decoh_times.show_matrix()

        print("Decoherence times matrix (fs):")
        tmp = decoh_times * units.au2fs
        tmp.show_matrix()

        print("Decoherence rates matrix (a.u.^-1):")
        decoh_rates.show_matrix()

    return decoh_times, decoh_rates


def decoherence_times_ave_old(Hvib, itimes, nsteps, verbosity=0):
    """

    Compute the matrix of decoherence times from the time-series data 
    that consists of vereral data sets

    Ref: Akimov, A. V; Prezhdo O. V. J. Phys. Chem. Lett. 2013, 4, 3857  


    Args:
        Hvib ( list of lists of CMATRIX objects ):
            Vibronic Hamiltonians along the trajectory for different data sets (adiabatic MDs),
            and potential different sections of the time-range
            where Hvib[idata][istep] is a CMATRIX object that represents a vibronic Hamiltonian
            from the data set ```idata``` at the time step ```istep```. 

        itimes ( list if ints ): initial times for averaging. The number ```itimes[idata]``` tells 
            which datapoint (timestep) of the time-series Hvib[idata] consider the beginning of 
            the range that will be used to compute gap fluctuations

        nsteps ( int ): the length of the time-steps to consider in the calculations for each data set

    Returns:
        tuple:  ( decoh_times, decoh_rates ), where

            * decoh_times ( MATRIX(N,N) ): the matrix of decoherence times, diagonal elements are
                set to very large number which corresponds to no decoherence of a state with itself.
                This value is computed based on gaps averaged over several data sets. [ units: a.u. ]

            * decoh_rates ( MATRIX(N,N) ): the matrix of decoherence rates, diagonal elements are
                set to zero which corresponds to no decoherence of a state with itself, the off-diagonal
                elements are equal to inverse of the off-diagonal matrix elements of ```decoh_times```.
                This value is computed based on gaps averaged over several data sets. [ units: a.u.^-1 ]

    """


    # Compute energy gaps
    dE = energy_gaps_ave(Hvib, itimes, nsteps)
    dE_ave, dE_std, dE_dw_bound, dE_up_bound = data_stat.mat_stat(dE)

    nstates = Hvib[0][0].num_of_cols
    decoh_times = MATRIX(nstates, nstates)   

    for a in range(0,nstates):
        for b in range(0,nstates):
            if a==b:
                decoh_times.set(a,a, 1.0e+10)
            else:
                de = dE_std.get(a,b)
                if de>0.0:
                      tau = math.sqrt(12.0/5.0) / de
                      decoh_times.set(a,b, tau)

    decoh_rates = decoherence_times2rates(decoh_times)

    if verbosity>0:
        print("Decoherence times matrix (a.u. of time):")
        decoh_times.show_matrix()

        print("Decoherence times matrix (fs):")
        tmp = decoh_times * units.au2fs
        tmp.show_matrix()

        print("Decoherence rates matrix (a.u.^-1):")
        decoh_rates.show_matrix()

    return decoh_times, decoh_rates



def decoherence_times_ave(Hvib, itimes, nsteps, verbosity=0):
    """

    Compute the matrix of decoherence times from the time-series data 
    that consists of vereral data sets

    Ref: Akimov, A. V; Prezhdo O. V. J. Phys. Chem. Lett. 2013, 4, 3857  


    Args:
        Hvib ( list of lists of CMATRIX objects ):
            Vibronic Hamiltonians along the trajectory for different data sets (adiabatic MDs),
            and potential different sections of the time-range
            where Hvib[idata][istep] is a CMATRIX object that represents a vibronic Hamiltonian
            from the data set ```idata``` at the time step ```istep```. 

        itimes ( list if ints ): initial times for averaging. The number ```itimes[idata]``` tells 
            which datapoint (timestep) of the time-series Hvib[idata] consider the beginning of 
            the range that will be used to compute gap fluctuations

        nsteps ( int ): the length of the time-steps to consider in the calculations for each data set

    Returns:
        tuple:  ( decoh_times, decoh_rates ), where

            * decoh_times ( MATRIX(N,N) ): the matrix of decoherence times, diagonal elements are
                set to very large number which corresponds to no decoherence of a state with itself.
                This value is computed based on gaps averaged over several data sets. [ units: a.u. ]

            * decoh_rates ( MATRIX(N,N) ): the matrix of decoherence rates, diagonal elements are
                set to zero which corresponds to no decoherence of a state with itself, the off-diagonal
                elements are equal to inverse of the off-diagonal matrix elements of ```decoh_times```.
                This value is computed based on gaps averaged over several data sets. [ units: a.u.^-1 ]

    """


    # Compute energy gaps
    ndata = len(Hvib)
    nitimes = len(itimes)
    nstates = Hvib[0][0].num_of_cols

    # Compute a concatenated list of gap magnitudes for all sub-trajectories    
    dE = []
    for idata in range(0,ndata):
        for it_indx in range(0,nitimes):
            it = itimes[it_indx]

            #========== Fluctuations for a given sub-trajectory ===========
            de = []
            for step in range(0,nsteps):
                dEij = MATRIX(nstates, nstates)

                for i in range(0,nstates):
                    for j in range(i+1, nstates):

                        deij = math.fabs(Hvib[idata][it+step].get(i,i).real - Hvib[idata][it+step].get(j,j).real)
                        dEij.set(i,j, deij)
                        dEij.set(j,i, deij)

                de.append(dEij)

            #============== Averages for this sub-trajectory ==================
            # we only care about dE_ave, other outputs are not used here
            dE_ave, dE_std, dE_dw_bound, dE_up_bound = data_stat.mat_stat(de)

            #============== Subtract the sub-trajectory-specific average =============
            for step in range(0,nsteps):
                dE.append( de[step] - dE_ave)


    # Compute the statistics of the obtained data set
    # the dE_ave should be zero, considering the transformations above
    # we only care about dE_std, other outputs are not used here
    dE_ave, dE_std, dE_dw_bound, dE_up_bound = data_stat.mat_stat(dE)


    decoh_times = MATRIX(nstates, nstates)   

    for a in range(0,nstates):
        for b in range(0, nstates):
            if a==b:
                decoh_times.set(a,a, 1.0e+10)
            else:
                de = dE_std.get(a,b)
                if de>0.0:
                      tau = math.sqrt(12.0/5.0) / de
                      decoh_times.set(a,b, tau)

    decoh_rates = decoherence_times2rates(decoh_times)

    if verbosity>0:
        print("Decoherence times matrix (a.u. of time):")
        decoh_times.show_matrix()

        print("Decoherence times matrix (fs):")
        tmp = decoh_times * units.au2fs
        tmp.show_matrix()

        print("Decoherence rates matrix (a.u.^-1):")
        decoh_rates.show_matrix()

    return decoh_times, decoh_rates



