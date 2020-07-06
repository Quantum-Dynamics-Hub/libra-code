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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import libra_py.common_utils as comn
import util.libutil as comn

import libra_py.data_read as data_read
from . import decoherence_times as dectim
import libra_py.tsh as tsh
import libra_py.tsh_stat as tsh_stat
import libra_py.units as units


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

    Hvib = []
    for i in range(0,params["nfiles"]):

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


