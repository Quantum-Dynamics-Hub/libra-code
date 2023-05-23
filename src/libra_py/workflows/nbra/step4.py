#*********************************************************************************
#* Copyright (C) 2017-2021 Brendan A. Smith, Wei Li, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
#
#  
#
"""

.. module:: step4
   :platform: Unix, Windows
   :synopsis: This module implements functions for running TSH calculations within the NBRA
       Example:    

       # Regular NBRA-NA-MD:

       >>> import step4
       >>> Hvib = step4.get_Hvib2(params)  # get the Hvib for all data sets
       >>> step4.transform_data(Hvib, {})  # default parameters don't change data
       >>> step4.run(Hvib, params)         # this ```params``` could be the same or different from the above

       # On-the-fly QSH-NA-MD:

       >>> import step4
       >>> import qsh
       >>> Hvib = qsh.run(qsh_params)      # generate the QSH Hvib data sets, see the ```qsh.run()``` for 
       >>>                                 # the description of the ```qsh_params``` parameters
       >>> step4.transform_data(Hvib, {})  # default parameters don't change data
       >>> step4.run(Hvib, params)         # 


.. moduleauthor:: Brendan A. Smith, Wei Li, Alexey V. Akimov


"""

import sys
import cmath
import math
import os
import multiprocessing as mp
import time
import numpy as np
import scipy.sparse as sp

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import libra_py.common_utils as comn
import util.libutil as comn

import libra_py.data_read as data_read
from . import decoherence_times as dectim
#import libra_py.tsh as tsh
import libra_py.tsh_stat as tsh_stat
import libra_py.units as units
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.recipes as recipes
import libra_py.fit as fit


def get_Hvib(params):
    """Read a single set of vibronic Hamiltonian files 

    Args:
        params ( dictionary ): parameters controlling the function execution

            Required parameter keys:

            * **params["nstates"]** ( int ): how many lines/columns in the file [Required!]
            * **params["nfiles"]** ( int ): how many files to read, starting from index 0 [Required!]
            * **params["Hvib_re_prefix"]** ( string ): prefixes of the files with real part of the Hvib(t) [Required!]
            * **params["Hvib_im_prefix"]** ( string ): prefixes of the files with imaginary part of the Hvib(t) [Required!]
            * **params["active_space"]** ( list of ints ): the indices of the states we care 
                about. These indices will be used to determine the size of the created CMATRIX objects
                and only these states will be extracted from the original files [ default: range(nstates) ]
            * **params["Hvib_re_suffix"]** ( string ): suffixes of the files with real part of the Hvib(t) [default: "_re"]
            * **params["Hvib_im_suffix"]** ( string ): suffixes of the files with imaginary part of the Hvib(t) [default: "_im"]

    Returns:
        list of CMATRIX objects: Hvib: 
            a time series of Hvib matrices, such that Hvib[time] is a Hvib at time step `time`

    Example:
        This example will read 10 pairs of files: "Hvib_0_re", "Hvib_0_im", "Hvib_1_re", "Hvib_1_im", ...
        "Hvib_9_re", "Hvib_9_im". Each file should contain a 4 x 4 matrix of numbers. It will generate a 
        list of 4 x 4 complex-valued matrices.

        >>> hvib = get_Hvib({"nstates":4, "nfiles":10, "Hvib_re_prefix":"Hvib", "Hvib_im_prefix":"Hvib"})


        The following example will do the same as the example above, however the intially-read 4 x 4 matrices will
        be partially discarded. Out of 16 values only 4 (the upper left block of 4 numbers)  will be stored in 
        the resulting list of 2 x 2 complex-valued matrices. 

        >>> hvib = get_Hvib({"nstates":4, "nfiles":10, "Hvib_re_prefix":"Hvib", "Hvib_im_prefix":"Hvib", "active_space":[0,1]})


    """

    critical_params = ["nstates", "nfiles", "Hvib_re_prefix", "Hvib_im_prefix"]
    default_params = { "Hvib_re_suffix":"_re", "Hvib_im_suffix":"_im", "active_space":range(params["nstates"])}
    comn.check_input(params, default_params, critical_params)

    nstates = params["nstates"]  # the number of states in the input files

    # init_time
    init_time = params["init_times"][0]

    Hvib = []
    # The previous was based on nfiles now it is based on the initial time and starts from init_time
    # for i in range(0,params["nfiles"]):
    for i in range(init_time, init_time+params["nfiles"]):

        filename_re = params["Hvib_re_prefix"]+str(i)+params["Hvib_re_suffix"]
        filename_im = params["Hvib_im_prefix"]+str(i)+params["Hvib_im_suffix"]
        hvib = data_read.get_matrix(nstates, nstates, filename_re, filename_im, params["active_space"] ) 
        Hvib.append(hvib)

    return Hvib


def get_Hvib2(params):
    """Reads several sets of vibronic Hamiltonian files 

    Args:
        params ( dictionary ): parameters controlling the function execution [Required!]

            Required parameter keys:

            * **params["data_set_paths"]** ( list of strings ):
                define the paths of the directories where the vibronic Hamiltonian files for
                different data sets (e.g. independent MD trajectories) are located. 
            .. note::
                In addition, requires parameters described in
                :func:`libra_py.workflows.nbra.step4.getHvib`

    Returns:
        list of lists of CMATRIX: Hvib: 
            the time series of Hvib matrices for several data sets, such that
            Hvib[idata][time] is a CMATRIX for the data set indexed by `idata`
            at time `time`

    Example:
        The full name of the vibronic Hamiltonian files read by this module should be:
    
        params["data_set_paths"][idata]+params["Hvib_re_prefix"]+integer(time step)+params["Hvib_re_suffix"] - for real part

        params["data_set_paths"][idata]+params["Hvib_im_prefix"]+integer(time step)+params["Hvib_im_suffix"] - for imaginary part

        Say, the directory "/home/alexeyak/test/step3/res0" contains files:
        Hvib_0_re, Hvib_1_re, .... ,    Hvib_999_re
        Hvib_0_im, Hvib_1_im, .... ,    Hvib_999_im

        Then set:

        >>> params["data_set_paths"] = ["/home/alexeyak/test/step3/res0/"]
        >>> params["Hvib_re_prefix"] = "Hvib_"
        >>> params["Hvib_re_suffix"] = "_re"
        >>> params["Hvib_im_prefix"] = "Hvib_"
        >>> params["Hvib_im_suffix"] = "_im"

    """

    critical_params = [ "data_set_paths" ] 
    default_params = { }
    comn.check_input(params, default_params, critical_params)

    H_vib = []

    for idata in params["data_set_paths"]:   # over all MD trajectories (data sets)
        prms = dict(params)    
        prms.update({"Hvib_re_prefix": idata+params["Hvib_re_prefix"] })
        prms.update({"Hvib_im_prefix": idata+params["Hvib_im_prefix"] })                

        h_vib = get_Hvib(prms)  
        H_vib.append(h_vib)

    return H_vib



def traj_statistics(i, Coeff, istate, Hvib, itimes):
    """Compute the averages over the TSH-ensembles

    Args:
        i ( int ): timestep index, counting since the beginning of the current sub-trajectory
        Coeff ( list of ntraj CMATRIX(nstates,1) object ): the TD-SE amplitudes for all trajectories
            (data sets/initial times/stochastic realizations)

        istate ( list of ntraj integers ): indices of the active states for each trajectory 
            (data sets/initial times/stochastic realizations)

        Hvib ( list of lists of CMATRIX(nstates,nstates) ): Hamiltonians for all data sets 
            and all (not just a sub-set of data!) timesteps
        itimes ( list of ints ): indices of the NA-MD starting points (in the global data indexing scale)

    Returns: 
        MATRIX(1, 3*nstates+4): the trajectory (and initial-condition)-averaged observables,
            the assumed format is: 

            First state info        ...          N-st state info All-states-related data

            E(0), P_SE(0), P_SH(0), ...,   E(nst-1), P_SE(nst-1), P_SH(nst-1), <E*P_SE>, <E*P_SH>, sum{P_SE}, sum{P_SH}
    
    """

    #================ Dimensions ==================
    Ntraj = len(Coeff)               # total number of trajectory = data set size x number of init times x number of SH trajectories
    nstates = Coeff[0].num_of_rows   # the number of states
    ndata = len(Hvib)                # how many data sets
    nitimes = len(itimes)            # how many initial times
    ntraj = int(Ntraj/(ndata * nitimes))  # how many stochastic SH trajectories per data set/initial condition


    # Update SE-derived density matrices
    denmat_se = tsh_stat.amplitudes2denmat(Coeff)   # list of Ntraj CMATRIX elements 

    # Use SE-derived density matrices to make SH-derived density matrices
    denmat_sh = []
    H_vib = []
    H_vib_ave = CMATRIX(nstates,nstates)  # Hvib averaged over the data sets/initial times

    for idata in range(0,ndata):
        for it_indx in range(0,nitimes):
            it = itimes[it_indx]

            H_vib_ave = H_vib_ave + Hvib[idata][it+i]

            for tr in range(0,ntraj):                
                Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                denmat_sh.append(CMATRIX(nstates, nstates))
                denmat_sh[Tr].set(istate[Tr],istate[Tr], 1.0, 0.0)
                H_vib.append(CMATRIX(Hvib[idata][it+i]))   

    H_vib_ave *= (1.0/float(ndata * nitimes))

    # Update TSH-ensemble-averaged SE and SH populations 
    ave_pop_sh, ave_pop_se = tsh_stat.ave_pop(denmat_sh, denmat_se)
    ave_en_sh,  ave_en_se  = tsh_stat.ave_en(denmat_sh, denmat_se, H_vib)

    # Save the computed data into a matrix to be output
    res = MATRIX(1, 3*nstates+4) 
   
    tot_sh, tot_se = 0.0, 0.0
    for j in range(0,nstates):
        res.set(0, 3*j+0, H_vib_ave.get(j,j).real)   # Energy of the state j
        res.set(0, 3*j+1, ave_pop_se.get(j,j).real)  # SE population
        res.set(0, 3*j+2, ave_pop_sh.get(j,j).real)  # SH population

        tot_se += ave_pop_se.get(j,j).real
        tot_sh += ave_pop_sh.get(j,j).real

    res.set(0, 3*nstates+0, ave_en_se)  # Average SE energy
    res.set(0, 3*nstates+1, ave_en_sh)  # Average SH energy
    res.set(0, 3*nstates+2, tot_se)     # Total SE population
    res.set(0, 3*nstates+3, tot_sh)     # Total SH population

    return res




def traj_statistics2(i, Pop, istate, Hvib, itimes):
    """Compute the averages over the TSH-ensembles

    Args:
        i ( int ): timestep index, counting since the beginning of the current sub-trajectory
        Pop ( list of ntraj CMATRIX(nstates, 1) object ): the qunatum populations for all trajectories
            (data sets/initial times/stochastic realizations)

        istate ( list of ntraj integers ): indices of the active states for each trajectory 
            (data sets/initial times/stochastic realizations)

        Hvib ( list of lists of CMATRIX(nstates,nstates) ): Hamiltonians for all data sets 
            and all (not just a sub-set of data!) timesteps
        itimes ( list of ints ): indices of the NA-MD starting points (in the global data indexing scale)

    Returns: 
        MATRIX(1, 3*nstates+4): the trajectory (and initial-condition)-averaged observables,
            the assumed format is: 

            First state info        ...          N-st state info All-states-related data

            E(0), P_SE(0), P_SH(0), ...,   E(nst-1), P_SE(nst-1), P_SH(nst-1), <E*P_SE>, <E*P_SH>, sum{P_SE}, sum{P_SH}
    
    """

    #================ Dimensions ==================
    Ntraj = len(Pop)                 # total number of trajectory = data set size x number of init times x number of SH trajectories
    nstates = Pop[0].num_of_rows     # the number of states
    ndata = len(Hvib)                # how many data sets
    nitimes = len(itimes)            # how many initial times
    ntraj = int(Ntraj/(ndata * nitimes))  # how many stochastic SH trajectories per data set/initial condition


    # Update SE-derived density matrices
    denmat_se = tsh_stat.pops2denmat(Pop)   # list of Ntraj CMATRIX(nstates,nstates) elements 

    # Update the SH-derived density matrices
    denmat_sh = []
    H_vib = []
    H_vib_ave = CMATRIX(nstates,nstates)  # Hvib averaged over the data sets/initial times

    for idata in range(0,ndata):
        for it_indx in range(0,nitimes):
            it = itimes[it_indx]

            H_vib_ave = H_vib_ave + Hvib[idata][it+i]

            for tr in range(0,ntraj):                
                Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                denmat_sh.append(CMATRIX(nstates, nstates))
                denmat_sh[Tr].set(istate[Tr],istate[Tr], 1.0, 0.0)
                H_vib.append(CMATRIX(Hvib[idata][it+i]))   

    H_vib_ave *= (1.0/float(ndata * nitimes))

    # Update TSH-ensemble-averaged SE and SH populations 
    ave_pop_sh, ave_pop_se = tsh_stat.ave_pop(denmat_sh, denmat_se)
    ave_en_sh,  ave_en_se  = tsh_stat.ave_en(denmat_sh, denmat_se, H_vib)

    # Save the computed data into a matrix to be output
    res = MATRIX(1, 3*nstates+4) 
   
    tot_sh, tot_se = 0.0, 0.0
    for j in range(0,nstates):
        res.set(0, 3*j+0, H_vib_ave.get(j,j).real)   # Energy of the state j
        res.set(0, 3*j+1, ave_pop_se.get(j,j).real)  # SE population
        res.set(0, 3*j+2, ave_pop_sh.get(j,j).real)  # SH population

        tot_se += ave_pop_se.get(j,j).real
        tot_sh += ave_pop_sh.get(j,j).real

    res.set(0, 3*nstates+0, ave_en_se)  # Average SE energy
    res.set(0, 3*nstates+1, ave_en_sh)  # Average SH energy
    res.set(0, 3*nstates+2, tot_se)     # Total SE population
    res.set(0, 3*nstates+3, tot_sh)     # Total SH population

    return res




def traj_statistics2_fast(i, Pop, istate, Hvib, itimes):
    """Compute the averages over the TSH-ensembles

    This version is optimized by removing all the dancing around data 
    and doing the direct calculations asap - no intermediate massive memory allocation!

    Args:
        i ( int ): timestep index, counting since the beginning of the current sub-trajectory
        Pop ( list of ntraj CMATRIX(nstates, 1) object ): the quantum populations for all trajectories
            (data sets/initial times/stochastic realizations)

        istate ( list of ntraj integers ): indices of the active states for each trajectory 
            (data sets/initial times/stochastic realizations)

        Hvib ( list of lists of CMATRIX(nstates,nstates) ): Hamiltonians for all data sets 
            and all (not just a sub-set of data!) timesteps
        itimes ( list of ints ): indices of the NA-MD starting points (in the global data indexing scale)

    Returns: 
        MATRIX(1, 3*nstates+4): the trajectory (and initial-condition)-averaged observables,
            the assumed format is: 

            First state info        ...          N-st state info All-states-related data

            E(0), P_SE(0), P_SH(0), ...,   E(nst-1), P_SE(nst-1), P_SH(nst-1), <E*P_SE>, <E*P_SH>, sum{P_SE}, sum{P_SH}
    
    """

    #================ Dimensions ==================
    Ntraj = len(Pop)                 # total number of trajectory = data set size x number of init times x number of SH trajectories
    nstates = Pop[0].num_of_rows     # the number of states
    ndata = len(Hvib)                # how many data sets
    nitimes = len(itimes)            # how many initial times
    ntraj = int(Ntraj/(ndata * nitimes))  # how many stochastic SH trajectories per data set/initial condition


    E_adi_ave =  MATRIX(nstates, 1) # E for each state averaged over the data sets/initial times
    ave_pop_sh = MATRIX(nstates, 1) # SH populations averaged over the data sets/initial times
    ave_pop_se = MATRIX(nstates, 1) # SE populations averaged over the data sets/initial times
    ave_en_sh = 0.0  # SH-population-weighted energy
    ave_en_se = 0.0  # SE-population-weighted energy


    # Temporary
    ene = MATRIX(nstates, 1)
    pop = MATRIX(nstates, 1)

    for idata in range(0,ndata):
        for it_indx in range(0,nitimes):
            it = itimes[it_indx]

            for st in range(nstates):
                ene.set(st, 0, Hvib[idata][it+i].get(st, st).real )

            E_adi_ave += ene

            for tr in range(0,ntraj):                
                Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                ave_pop_sh.add(istate[Tr], 0, 1.0) 
                ave_en_sh += ene.get(istate[Tr], 0)

                pop = Pop[Tr].real()
                ave_pop_se += pop
                ave_en_se += (ene.T() * pop).get(0,0)
                

    nrm1 = (1.0/float(ndata * nitimes))
    nrm2 = nrm1/float(ntraj)

    E_adi_ave *= nrm1
    ave_pop_sh *= nrm2
    ave_pop_se *= nrm2
    ave_en_sh *= nrm2
    ave_en_se *= nrm2


    # Save the computed data into a matrix to be output
    res = MATRIX(1, 3*nstates+4) 
   
    tot_sh, tot_se = 0.0, 0.0
    for j in range(0,nstates):
        res.set(0, 3*j+0, E_adi_ave.get(j,0))   # Energy of the state j
        res.set(0, 3*j+1, ave_pop_se.get(j,0))  # SE population
        res.set(0, 3*j+2, ave_pop_sh.get(j,0))  # SH population

        tot_se += ave_pop_se.get(j,0)
        tot_sh += ave_pop_sh.get(j,0)

    res.set(0, 3*nstates+0, ave_en_se)  # Average SE energy
    res.set(0, 3*nstates+1, ave_en_sh)  # Average SH energy
    res.set(0, 3*nstates+2, tot_se)     # Total SE population
    res.set(0, 3*nstates+3, tot_sh)     # Total SH population

    return res




def printout(t, res, outfile):
    """This function does a simple output of a matrix columns to a file

    Args:    
        t ( double ): time [units: a.u.] 
        res ( MATRIX(1,N) ): information to be printed out
        outfile ( string ): filename where we'll print everything out, the output 
            will be appended to the existing output file

    Returns:
        None: but modifies the file

    """

    N = res.num_of_cols

    line = "%8.5f " % (t)
    for i in range(0,N):
        line = line + " %8.5f " % (res.get(0,i))
    line = line + " \n"

    f = open(outfile, "a") 
    f.write(line)
    f.close()



# This function is deprecated and does not work in the current version of Libra, 5/18/2023
def run(H_vib, params):
    """
    
    The main procedure to run NA-MD calculations within the NBRA workflow

    Args: 
        H_vib ( list of lists of CMATRIX objects ): the vibronic Hamiltonian for all data sets and all time-points
            H_vib[idata][istep].get(i,j) - i,j matrix element for the data set `idata` and step in that data set `istep`
    
        params ( dictionary ): the parameters that control the execution of the NA-MD-NBRA calculations

            * **params["nsteps"]** ( int ): the length of the NA-MD trajectory. This parameter is not 
                necessarily the same as len(H_vib[0]), so need to be provided [Required!] 

            * **params["T"]** ( double ): temperature of nuclear/electronic dynamics [in K, default: 300.0]
            * **params["ntraj"]** ( int ): the number of stochastic surface hopping trajectories [default: 1]
            * **params["tdse_Ham"]** ( int ): option to select either the regular (input) or Boltzmann-corrected
                Hamiltonian:

                - 0 - regular [ default ]
                - 1 - Boltzmann-corrected 

            * **params["Hvib_type"]** ( int ): option to select if the Hvib is a diabatic or an adiabatic
                Hamiltonian:

                - 0 - diabatic
                - 1 - adiabatic  [ default ]

            * **params["sh_method"]** ( int ): selects the algorithm to compute surface hopping probabilities 
                Options:

                - 0 - MSSH
                - 1 - FSSH [ default ]

            * **params["decoherence_constants"]** ( int ): selects whether to compute decoherence parameters
                on the fly or to use provided parameters:

                - 0 - pre-compute the parameters from the trajectory data before NA-MD run [ default ]
                - 1 - use the provided parameters ..seealso:: ```params["decoherence_times"]``` and ```params["decoherence_rates"]```
                - 20 - use the time-dependent decoherence times as in DISH paper. This is different from
                    option 0 in that these numbers depend on the state amplitudes. Dephasing times are computed as in 0.
                - 21 - use the time-dependent decoherence times as in DISH paper. This is different from
                    option 0 in that these numbers depend on the state amplitudes. Dephasing times are computed as in 1.


            * **params["decoherence_times"]** ( MATRIX(nstates,nstates) ): decoherence times for all 
                pairs of states. This should be provided if ``` params["decoherence_constants"] == 1``` the dimensions should be
                consistent with those of the input Hvib data. [ units: a.u. of time ]
             

            * **params["decoherence_method"]** ( int ): selects the decoherence method 
                Options:

                - 0 - no decoherence [ default ]
                - 1 - ID-A
                - 2 - MSDM
                - 3 - DISH

            * **params["dt"]** ( double ): nuclear dynamics integration time step [in a.u. of time, default: 41.0]
            * **params["Boltz_opt"]** ( int ): option to select a probability of hopping acceptance [default: 3]
                Options:

                - 0 - all proposed hops are accepted - no rejection based on energies
                - 1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
                - 2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
                - 3 - generalization of "1", but actually it should be changed in case there are many degenerate levels

            * **params["istate"]** ( int ): index of the initial state [default: 0]
            * **params["init_times"]** ( list of ints ): indices of the starting point in the provided data arrays [default: [0]]
            * **params["outfile"]** ( string ): the name of the file where to print populations
                and energies of states [default: "_out.txt"]    

    Returns: 
        MATRIX(nsteps, 3*nstates+5): the trajectory (and initial-condition)-averaged observables for every timesteps,
            the assumed format is: 

            time, first state info        ...          N-st state info All-states-related data

            time, E(0), P_SE(0), P_SH(0), ...,   E(nst-1), P_SE(nst-1), P_SH(nst-1), <E*P_SE>, <E*P_SH>, sum{P_SE}, sum{P_SH}

     
    """

    
    critical_params = [ "nsteps" ] 
    default_params = { "T":300.0, "ntraj":1,
                       "tdse_Ham":0, "sh_method":1, "decoherence_constants": 0, "decoherence_method":0, "dt":41.0, "Boltz_opt":3,
                       "Hvib_type":1,
                       "istate":0, "init_times":[0], "outfile":"_out.txt" }
    comn.check_input(params, default_params, critical_params)


    rnd = Random()

    ndata = len(H_vib)
    nsteps = params["nsteps"]
    nstates = H_vib[0][0].num_of_cols  # number of states

    ntraj = params["ntraj"]
    nitimes = len(params["init_times"])
    Ntraj = ndata * nitimes * ntraj

    T = params["T"]
    bolt_opt = params["Boltz_opt"]
    dt = params["dt"]
    tdse_Ham = params["tdse_Ham"]

    res = MATRIX(nsteps, 3*nstates+5)


    #========== Compute PARAMETERS  ===============
    # Decoherence times
    # these are actually the dephasing rates!
    tau, dephasing_rates = None, None

    if params["decoherence_constants"] == 0 or params["decoherence_constants"]==20:
        tau, dephasing_rates = dectim.decoherence_times_ave(H_vib, params["init_times"], nsteps, 1) 

    elif params["decoherence_constants"] == 1 or params["decoherence_constants"]==21:
        if params["decoherence_times"].num_of_cols != nstates:
            print("Error: dimensions of the input decoherence times matrix are not consistent with \
                   the dimensions of the Hamiltonian matrices (the number of states). Exiting...\n")
            sys.exit(0)
        else:
            tau = MATRIX(params["decoherence_times"])
            dephasing_rates = dectim.decoherence_times2rates(tau)



    #========== Initialize the DYNAMICAL VARIABLES  ===============
    # TD-SE coefficients and active state indices
    Coeff, istate = [], []

    # Coherence times and coherence intervals for DISH
    t_m, tau_m = [], []

    for tr in range(0,Ntraj):
        istate.append(params["istate"])
        Coeff.append(CMATRIX(nstates, 1)); 
        Coeff[tr].set(params["istate"], 1.0, 0.0)
        t_m.append(MATRIX(nstates,1))
        tau_m.append(MATRIX(nstates,1))  

           
    # Prepare the output file
    f = open(params["outfile"],"w"); f.close()

    #=============== Entering the DYNAMICS ========================
    for i in range(0,nsteps):  # over all evolution times


        #============== Analysis of the Dynamics  =================
        # Compute the averages
        res_i = traj_statistics(i, Coeff, istate, H_vib, params["init_times"])

        # Print out into a file
        printout(i*dt, res_i, params["outfile"])

        # Update the overal results matrix
        res.set(i,0, i*dt)
        push_submatrix(res, res_i, Py2Cpp_int([i]), Py2Cpp_int( list(range(1,3*nstates+5)) ) )


        #=============== Propagation ==============================
        for idata in range(0,ndata):   # over all MD trajectories (data sets)

            for it_indx in range(0,nitimes): # over all initial times

                it = params["init_times"][it_indx]

                for tr in range(0,ntraj):  # over all stochastic trajectories

                    Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                    #============== Propagation: TD-SE and surface hopping ==========
                    # Coherent evolution amplitudes
                    Heff = None 

                    if tdse_Ham==0:
                        Heff = H_vib[idata][it+i]
                    elif tdse_Ham==1:
                        Heff = tsh.Boltz_corr_Ham(H_vib[idata][it+i], Coeff[Tr], params["T"], params["Hvib_type"])

                    propagate_electronic(dt, Coeff[Tr], Heff)   # propagate the electronic DOFs

        
                    # Surface hopping 
                    ksi  = rnd.uniform(0.0, 1.0)
                    ksi2 = rnd.uniform(0.0, 1.0)
        
                    if params["decoherence_method"] in [0, 1, 2]:
                    
                        do_collapse = None
                        if params["decoherence_method"]==0:    # No decoherence
                            do_collapse = 0
                        elif params["decoherence_method"]==1:  # ID-A, taken care of in the tsh.hopping
                            do_collapse = 1
                        elif params["decoherence_method"]==2:  # MSDM
                            do_collapse = 0                            
                            Coeff[Tr] = sdm(Coeff[Tr], dt, istate[Tr], dephasing_rates)
 
                        istate[Tr], Coeff[Tr] = tsh.hopping(Coeff[Tr], Heff, istate[Tr], params["sh_method"], do_collapse, ksi, ksi2, dt, T, bolt_opt)
                    
                    elif params["decoherence_method"] in [3]:  # DISH
                    
                        tau_m[Tr] = coherence_intervals(Coeff[Tr], dephasing_rates)
                        istate[Tr] = tsh.dish_py(Coeff[Tr], istate[Tr], t_m[Tr], tau_m[Tr], Heff, bolt_opt, T, ksi, ksi2)
                        t_m[Tr] += dt

        
    return res



def run_tsh(common_params, compute_model, model_params):    
    """
    This function sleeps for a while, initialized nuclear and electronic variables, 
    and executed the NA-MD calculations

    Args: 
        dyn_params ( dict ): parameters controlling the NAMD calculations as used by the `run_tsh` function above
            Note that this is only a template - some parameters will be modified internally according to the 
            other parameters of this function

        compute_model (Python function object): the function that computes the Hamiltonian and other properties

        model_params ( dict ): parameters of the `compute_model` function

    Returns:
        None: but produces file outputs

    """

    params = dict(common_params) 

    critical_params = [ ] 
    default_params = { "wait_time":0.0, "ntraj":1, "nstates":1, "istate":[1, 0],
                       "x0":[0.0], "p0":[0.0], "masses":[1.0], "k":[0.01], 
                       "nucl_init_type":2, "elec_init_type":3
                     }
    comn.check_input(params, default_params, critical_params)
  


    time.sleep( int(params["wait_time"]) )

    # Random numbers generator object
    rnd = Random()
    
    #============ Initialize dynamical variables ==================
    x0, p0, masses, k0 = params["x0"], params["p0"], params["masses"], params["k"]
    ntraj, nstates = params["ntraj"], params["nstates"]
    inucl = params["nucl_init_type"]
    ielec = params["elec_init_type"]
    
    # Nuclear
    init_nucl = {"init_type":inucl, "force_constant":k0, "ntraj":ntraj}
    q, p, iM = tsh_dynamics.init_nuclear_dyn_var(x0, p0, masses, init_nucl, rnd)
    
    # Electronic
    istate = params["istate"]
    istates = []
    for i in range(nstates):
        istates.append(0.0)
    istates[ istate[1] ] = 1.0    
    _init_elec = { "init_type":ielec, "nstates":nstates, "istates":istates, "rep":istate[0],  "ntraj":ntraj   }
    

    #============= Dynamical variables ==============
    dyn_params = dict(common_params)
    
    # This should update only the properties that aren't defined, but not override the existing values!
    critical_params = [  ]     
    default_params = { "mem_output_level":3 }     
    comn.check_input(dyn_params, default_params, critical_params)
                    
    _model_params = dict(model_params)
    _model_params.update({"model0": model_params["model"] })
               
    start = time.time()        
    res = tsh_dynamics.generic_recipe(q, p, iM, dyn_params, compute_model, _model_params, _init_elec, rnd)
    end = time.time()    
    print(F"Calculation time = {end - start} seconds")

    
def make_var_pool(dyn_params, compute_model, model_params, _rnd, 
                  method_names_map = {0:"FSSH", 1:"IDA", 2:"mSDM", 3:"DISH", 21:"mSDM2", 31:"DISH2" },
                  init_states = [1], tsh_methods = [0], batches = list(range(10)), ham_rep=1, is_nbra=1 ):
    """
    This function prepares the variables pool for the parallelization 

    Args: 
        dyn_params ( dict ): parameters controlling the NAMD calculations as used by the `run_tsh` function above
            Note that this is only a template - some parameters will be modified internally according to the 
            other parameters of this function

        compute_model (Python function object): the function that computes the Hamiltonian and other properties

        model_params ( dict ): parameters of the `compute_model` function

        _rnd ( Random ): random numbers generator instance

        methods_names_map ( dict ): the mapping between the TSH method indices and their names - needed for 
             making the output directories [ default: {0:"FSSH", 1:"IDA", 2:"mSDM", 3:"DISH", 21:"mSDM2", 31:"DISH2" } ]

        init_states ( int list ): the indices of all initial states for which to do the computations [ default: [1] ]

        tsh_methods ( int list ): the indices of all methods for which to do the computations - must be consistent with 
             the `methods_names_map` variable [ default: [0] ]
 
        batches ( int list ): how many batches of simulations to run [ default: list(range(10)) ]

        ham_rep ( int ): what type of representation is used to compute the adiabatic Hamiltonian:

           - 0 - diabatic, it gets diagonalized to yield the adiabatic properties;
           - 1 - adiabatic properties are already known [default ]

        is_nbra ( int ): selector to choose whether we are doing an NBRA simulations or not
  
           - 0 - not an NBRA
           - 1 - NBRA [ default  ]

    Returns:
        list of tuples: each tuple contains the input variables for the function that we want run in multiple threads

    """
    
    # Create the pool of variables 
    var_pool = []
    for istate in init_states:  # initial states   
        for method in tsh_methods:  # decoherence method: FSSH, IDA, mSDM, DISH        
            
            name = method_names_map[method]   
            
            for batch in batches:        
                mdl_prms = dict(model_params)
                
                
                prms = dict(dyn_params)                
                                
                prms["wait_time"] = _rnd.uniform(0.0, 5.0)                                                
                dir_prefix = prms["dir_prefix"]
                                
                prms["prefix"] = F"{dir_prefix}/start_s{istate}_{name}_batch{batch}"                
                prms["prefix2"] = F"{dir_prefix}/_start_s{istate}_{name}_batch{batch}" 
                prms["istate"] = [1, istate] # adiabatic state `istate`; Recall index from 0

                # This instruction is made backward-compatible with the commented block
                recipes.set_method(prms, ham_rep, is_nbra, method)
                                
                var_pool.append( (prms, compute_model, mdl_prms) )    

    return var_pool



def namd_workflow(dyn_params, compute_model, model_params, _rnd, nthreads = 4,   
                  method_names_map = {0:"FSSH", 1:"IDA", 2:"mSDM", 3:"DISH", 21:"mSDM2", 31:"DISH2" },
                  init_states = [1], tsh_methods = [0], batches = list(range(10)), fork_or_spawn="fork", parallel=True,
                  ham_rep = 1, is_nbra = 1):
    """
    
    This is a new top-level wrapper to run NA-MD calculations within the NBRA workflow (although it could be more general too)
    This wrapper provides a multithreading parallelization over: initial state - TSH method - batches 

    Args: 
        dyn_params ( dict ): parameters controlling the NAMD calculations as used by the `run_tsh` function above
            Note that this is only a template - some parameters will be modified internally according to the 
            other parameters of this function

        compute_model (Python function object): the function that computes the Hamiltonian and other properties

        model_params ( dict ): parameters of the `compute_model` function

        _rnd ( Random ): random numbers generator instance

        nthreads ( int ): the number of threads over which to parallelize the calculations [ default: 4]

        methods_names_map ( dict ): the mapping between the TSH method indices and their names - needed for 
             making the output directories [ default: {0:"FSSH", 1:"IDA", 2:"mSDM", 3:"DISH", 21:"mSDM2", 31:"DISH2" } ]

        init_states ( int list ): the indices of all initial states for which to do the computations [ default: [1] ]

        tsh_methods ( int list ): the indices of all methods for which to do the computations - must be consistent with 
             the `methods_names_map` variable [ default: [0] ]
 
        batches ( int list ): how many batches of simulations to run [ default: list(range(10)) ]

        fork_or_spawn ( string ): what kind of multiprocessing to use. Can be of only 2 types:
           
           - fork : via forking [ default ]
           - spawn: via spawning

        parallel ( Boolean ): whether to run calculations with the multiprocessing parallelization [ default : True]

        ham_rep ( int ): what type of representation is used to compute the adiabatic Hamiltonian:

           - 0 - diabatic, it gets diagonalized to yield the adiabatic properties;
           - 1 - adiabatic properties are already known [default ]

        is_nbra ( int ): selector to choose whether we are doing an NBRA simulations or not
  
           - 0 - not an NBRA
           - 1 - NBRA [ default  ]


    Returns:
        None: but it call functions that produce outputs elsewhere
       
    """

               
    var_pool = make_var_pool(dyn_params, compute_model, model_params, _rnd, 
                             method_names_map, init_states, tsh_methods, batches, ham_rep, is_nbra )
    
    if parallel==True:    
        t1 = time.time()
        pool = mp.get_context(F"{fork_or_spawn}").Pool( nthreads )
    #    pool = mp.Pool( nthreads )
        pool.starmap( run_tsh, var_pool )
        pool.close()
        pool.join()

        t2 = time.time()
        print(F"Total time {t2 - t1}") 

    else:
        t1 = time.time()
        for var in var_pool:
            run_tsh(var[0], var[1], var[2])
        t2 = time.time()
        print(F"Total time {t2 - t1}") 

    

def nice_plots(dyn_params, init_states, tsh_methods, methods, batches, fig_label="NA-MD", txt_type=0):
    """
    This function produces nice plots of the ground state populations for all calculations
    """

    plt.rc('axes', titlesize=12)      # fontsize of the axes title
    plt.rc('axes', labelsize=12)      # fontsize of the x and y labels
    plt.rc('legend', fontsize=10)     # legend fontsize
    plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
    plt.rc('figure.subplot', left=0.2)
    plt.rc('figure.subplot', right=0.95)
    plt.rc('figure.subplot', bottom=0.13)
    plt.rc('figure.subplot', top=0.88)
    
    prefix = dyn_params["dir_prefix"]
    all_time = np.array( list(range( dyn_params["nsteps"] ) ) ) * dyn_params["dt"] * units.au2fs
    
    for istate in init_states:  # initial states       
        print(F"======= Running initial state {istate} =======")
        for method in tsh_methods:  # decoherence method: FSSH, IDA, mSDM, DISH                        
            name = methods[method]
            print(F"    *** Running decoherence method {method} ( {name} ) *** ")
            
            figure = plt.figure(num=None, figsize=(3.21, 2.41), dpi=300, edgecolor='black', frameon=True)        
            name = methods[method]        
                    
            tau, tau2 = [], []
            rat, rat2 = [], []
            for batch in batches:        
                            
                if txt_type==0:
                    infile1 = F"{prefix}/start_s{istate}_{name}_batch{batch}/time.txt"
                    infile2 = F"{prefix}/start_s{istate}_{name}_batch{batch}/SH_pop.txt"
                elif txt_type==1:
                    infile1 = F"{prefix}/_start_s{istate}_{name}_batch{batch}/time.txt"
                    infile2 = F"{prefix}/_start_s{istate}_{name}_batch{batch}/SH_pop.txt"

                                        
                # Plot the raw data
                t = np.array(data_read.get_data_from_file2(infile1, [0])[0]) * units.au2fs
                p = np.array(data_read.get_data_from_file2(infile2, [0,1]))
                
                print(F"raw data lengths: { len(t)} { len(p[0]) }")
                            
                if batch==0:
                    plt.plot(t, p[0], color="black", label="", linewidth=1)      
                    #plt.plot(t, p[1], color="red", label="P1", linewidth=1)      
                    #plt.plot(t, p[2], color="blue", label="P2", linewidth=1)      
                    #plt.plot(t, p[3], color="green", label="P3", linewidth=1)      
                else:
                    plt.plot(t, p[0], color="black", label="", linewidth=1)      
                    #plt.plot(t, p[0], color=colors[ clrs_index[batch] ], label="", linewidth=1)
                    
                    #plt.plot(t, p[1], color="red", label="", linewidth=1)      
                    #plt.plot(t, p[2], color="blue", label="", linewidth=1)      
                    #plt.plot(t, p[3], color="green", label="", linewidth=1)      
                
                            
                
                T = []
                P = []
                #for ia, a in enumerate(p[0]):
                #    T.append(t[ia])
                #    P.append(p[0, ia])
                    
                for ia, a in enumerate(p[0]):
                    if a > 0.01 and a < 0.99 and t[ia] > 0.0:
                        T.append(t[ia])
                        P.append(p[0, ia])
                                                    
                if sum(P) > 0.01:
                    
                    # Use a well-behaved subset of data for fitting                        
                    pop = 1.0 - np.array(P)
                
                    # Do the fitting and collect the results
                    verb, opt = 0, 0
                    
                    fit_res, fit_a, fit_b = None, None, None
                    if istate in [-1] and method in [0]:
                        fit_res, fit_a, fit_b = fit.fit_gau(T, pop, 0.0, verb, opt)  
                        r = fit_b                    
                        rat2.append( r )            # 1/tau^2
                        rat.append( math.sqrt(r) )  # 1/tau                                         
                        tau2.append( 1.0/r )           # tau^2
                        tau.append(1.0/math.sqrt(r) )  # tau                    
                        
                    else:
                        fit_res, fit_a, fit_b = fit.fit_exp(T, pop, 0.0, verb, opt)                  
                        r = fit_b                                        
                        rat2.append(r**2)            # 1/tau^2
                        rat.append( r )              # 1/tau                                        
                        tau2.append( 1.0/(r**2) )    # tau^2
                        tau.append(1.0/r )           # tau
                                                                                                    
                    print( "Fitting parameters : A = ", fit_a, " and B = ", fit_b, " 1/B = ", 1/fit_b)                                                
                else:
                    pass
            
            n = len(rat)
            
            if n > 0:
                print(F"Number of samples = {n}")
                
                # Prefactor for 95% confidence interval
                # The error bars reporting is :  average +/-  prefactor * standard_deviation
                prefactor = 1.96/math.sqrt(n)
                
                # Statistical analysis of the obtained timescales
                rat_stat = DATA(rat); rat_stat.Calculate_Estimators()
                rat2_stat = DATA(rat2); rat2_stat.Calculate_Estimators()            
                tau_stat = DATA(tau); tau_stat.Calculate_Estimators()
                tau2_stat = DATA(tau2); tau2_stat.Calculate_Estimators()            
                
                
                print(F"Lifetime from averaging individual lifetimes {tau_stat.ave} +/- { prefactor * tau_stat.sd }")
                print(F"Average rates  {rat_stat.ave} +/- { prefactor * rat_stat.sd }")  
                                                    
                
                # Plot the average fitting     
                Y, Y1, Y2 = None, None, None
                TAU, ERR = None, None
                if istate in [-1] and method in [0]:   
                    """
                    tau = <r2> ^ {-1/2}
                    d tau = - 1/2 <r2>^ {-3/2} * d( r2)
                    """
                    
                    TAU = 1.0/math.sqrt(rat2_stat.ave)
                    ERR = prefactor * 0.5 * TAU**3 * rat2_stat.sd # tau_stat.sd #( tau_stat.ave / rat_stat.ave) *rat_stat.sd
                                    
                    
                    print(F"Lifetime { TAU } +/- { ERR }")    
                    Y = 1.0 - np.exp(-(all_time**2 * (rat2_stat.ave) )  )
                    Y1 = 1.0 - np.exp(-(all_time**2 * (rat2_stat.ave - prefactor * rat2_stat.sd) )  )
                    Y2 = 1.0 - np.exp(-(all_time**2 * (rat2_stat.ave + prefactor * rat2_stat.sd) )  ) 
                    
                                                                
                else:   
                    """
                    tau = <r> ^ {-1}
                    d tau = - <r>^ {-2} * d(r)
                    """
                    
                    TAU = 1.0/rat_stat.ave
                    ERR = prefactor * TAU**2 * rat_stat.sd 
                    
                    print(F"Lifetime {TAU} +/- { ERR }")
                    Y = 1.0 - np.exp(-(all_time * rat_stat.ave ) ) 
                    Y1 = 1.0 - np.exp(-(all_time * (rat_stat.ave - prefactor * rat_stat.sd ) )  )
                    Y2 = 1.0 - np.exp(-(all_time * (rat_stat.ave + prefactor * rat_stat.sd) )  )                                
                                            
                        
                plt.plot(all_time , Y, color="red", 
                         label=F"{name}, {TAU/1000.0 : 3.1f} +/- {ERR/1000.0 : 3.1f} ps", 
                         linewidth=2 )
                plt.plot(all_time, Y1, color="red", label=F"", linewidth=1, linestyle='dashed' )
                plt.plot(all_time, Y2, color="red", label=F"", linewidth=1, linestyle='dashed' )
                
                                
                        
            plt.title(F"{fig_label}",fontsize=9.5)
            plt.legend(fontsize=6.75, ncol=1, loc='upper left')
            plt.xlabel('Time, fs',fontsize=10)
            plt.ylabel('Population',fontsize=10)
            plt.tight_layout()
            if txt_type==0:
                plt.savefig(F'{prefix}/start_s{istate}_{name}.png', dpi=300)
            elif txt_type==1:
                plt.savefig(F'{prefix}/_start_s{istate}_{name}.png', dpi=300)

            plt.show()                                                                                                            


def get_Hvib_scipy(params):
    """
    This function is used to get the set of Hvibs (or any other set of data like 
    St matrices) in scipy.sparse format and return them as CMATRIX format.
    Args:
        params (dictionary): parameters controlling the function execution [Required!]

            Required parameter keys:

            * **params["data_set_paths"]** ( list of strings ):
                define the paths of the directories where the vibronic Hamiltonian files for
                different data sets (e.g. independent MD trajectories) are located. 
            .. note::
                In addition, requires parameters described in
                :func:`libra_py.workflows.nbra.step4.getHvib`

    Returns:
        list of lists of CMATRIX: Hvib: 
            the time series of Hvib matrices for several data sets, such that
            Hvib[idata][time] is a CMATRIX for the data set indexed by `idata`
            at time `time`
    """
    hvib = []
    for i in range(len(params['data_set_paths'])):
        hvib.append([])
        for step in range(params['nfiles']):
            step += params['init_times']
    
            try:
                file_name = params['data_set_paths'][0] + params['Hvib_re_prefix'] + str(step) + params['Hvib_re_suffix']
                tmp_real = sp.load_npz(file_name).real
            except:
                print(F'File {file_name} not found! Please check the path to Hvibs. Setting zero matrix into hvib...')
                tmp_real = sp.csc_matrix( np.zeros((params['nstates'], params['nstates'])) )
    
            try:
                file_name = params['data_set_paths'][0] + params['Hvib_im_prefix'] + str(step) + params['Hvib_im_suffix']
                tmp_imag = sp.load_npz(file_name).real
            except:
                print(F'File {file_name} not found! Please check the path to Hvibs. Setting zero matrix into hvib...')
                tmp_imag = sp.csc_matrix( np.zeros((params['nstates'], params['nstates'])) )
    
            real_MATRIX = scipynpz2MATRIX(tmp_real)
            imag_MATRIX = scipynpz2MATRIX(tmp_imag)
            hvib[i].append( CMATRIX(real_MATRIX, imag_MATRIX)  )
    
    return hvib

   
