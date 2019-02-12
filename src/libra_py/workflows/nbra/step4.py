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
       >>> 
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

import common_utils as comn
import decoherence_times as dectim
import libra_py.tsh as tsh
import libra_py.units as units
import qsh


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
        hvib = comn.get_matrix(nstates, nstates, filename_re, filename_im, params["active_space"] ) 
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


def transform_data(X, params):
    """
    This is an auxiliary function to morph the original X (e.g. H_vib) according to:     

    X(original) ->    ( X + shift1 ) x (scale) + shift2,  

    Here, x indicates the element-wise multiplicaiton, and shift1, shift2, and scale are matrices

    Transformes X (e.g. H_vib) input directly, so changes the original input

    === Data tuning (e.g. units conversion or other trickery) ===
    params["shift1"]      [CMATRIX(nstates,nstates)] - first shift corrections [units of X]
    params["shift2"]      [CMATRIX(nstates,nstates)] - second shift corrections [units of X]
    params["scale"]       [CMATRIX(nstates,nstates)] - scaling of the X [unitless]


    """

    ndata = len(X)
    nsteps = len(X[0])
    nstates = X[0][0].num_of_cols

    sh1 = CMATRIX(nstates, nstates) # zero
    sh2 = CMATRIX(nstates, nstates) # zero
    scl = CMATRIX(nstates, nstates) # all elements are 1

    for i in xrange(nstates):
        for j in xrange(nstates):
            scl.set(i,j, 1.0+0.0j)

    critical_params = [  ] 
    default_params = { "shift1":sh1, "shift2":sh2, "scale":scl  }
    comn.check_input(params, default_params, critical_params)

    for idata in xrange(ndata):
        for istep in xrange(nsteps):

            tmp = CMATRIX(X[idata][istep])
            tmp = tmp + params["shift1"]
            tmp.dot_product( tmp, params["scale"] )
            tmp = tmp + params["shift2"]
            X[idata][istep] = CMATRIX(tmp)




def traj_statistics(i, Coeff, istate, Hvib, itimes):
    """
    Compute the averages over the TSH-ensembles

    i           [int]                         - timestep index
    Coeff       [list of CMATRIX(nstates,1)]  - TD-SE amplitudes for all trajectories (data sets/initial times/stochastic realizations)
    istate      [list of nstates integers]    - indices of the active states for each trajectory (data sets/initial times/stochastic realizations)
    Hvib        [list of lists of CMATRIX(nstates,nstates)]
                                              - Hamiltonians for all data sets and all timesteps
    itimes      [list of ints]                - indices of the NA-MD starting points


    Returns: 
    MATRIX(1, 3*nstates+4)   with the following columns: 
    E(0), P_SE(0), P_SH(0), ...,   E(nst-1), P_SE(nst-1), P_SH(nst-1), <E*P_SE>, <E*P_SH>, sum{P_SE}, sum{P_SH}
    
       First state info                    Nst state info

    """

    #================ Dimensions ==================

    #print "In traj_statistics...\n"

    Ntraj = len(Coeff)               # total number of trajectory = data set size x number of init times x number of SH trajectories
    nstates = Coeff[0].num_of_rows   # the number of states
    ndata = len(Hvib)                # how many data sets
    nitimes = len(itimes)            # how many initial times
    ntraj = Ntraj/(ndata * nitimes)  # how many stochastic SH trajectories per data set/initial condition

    #print "Ntraj = ", Ntraj
    #print "nstates = ", nstates
    #print "ndata = ", ndata
    #print "nitimes = ", nitimes
    #print "ntraj = ", ntraj


    # Update SE-derived density matrices
    denmat_se = tsh.amplitudes2denmat(Coeff)   # list of Ntraj CMATRIX elements 

    # Use SE-derived density matrices to make SH-derived density matrices
    denmat_sh = []
    H_vib = []
    H_vib_ave = CMATRIX(nstates,nstates)  # Hvib averaged over the data sets/initial times

    for idata in xrange(ndata):
        for it_indx in xrange(nitimes):
            it = itimes[it_indx]

            H_vib_ave = H_vib_ave + Hvib[idata][it+i]

            for tr in xrange(ntraj):                
                Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                denmat_sh.append(CMATRIX(nstates, nstates))
                denmat_sh[Tr].set(istate[Tr],istate[Tr], 1.0, 0.0)
                H_vib.append(CMATRIX(Hvib[idata][it+i]))   

    H_vib_ave *= (1.0/float(ndata * nitimes))

    # Update TSH-ensemble-averaged SE and SH populations 
    #print "len(denmat_se) = ", len(denmat_se)
    #print "len(denmat_sh) = ", len(denmat_sh)
    #print "len(H_vbi) = ", len(H_vib)

    ave_pop_sh, ave_pop_se = tsh.ave_pop(denmat_sh, denmat_se)
    ave_en_sh,  ave_en_se  = tsh.ave_en(denmat_sh, denmat_se, H_vib)

    # Save the computed data into a matrix to be output
    res = MATRIX(1, 3*nstates+4) 
   
    tot_sh, tot_se = 0.0, 0.0
    for j in xrange(nstates):
        res.set(0, 3*j+0, H_vib_ave.get(j,j).real)   # Energy of the state j
        res.set(0, 3*j+1, ave_pop_se.get(j,j).real)  # SE population
        res.set(0, 3*j+2, ave_pop_sh.get(j,j).real)  # SH population

        tot_se += ave_pop_se.get(j,j).real
        tot_sh += ave_pop_sh.get(j,j).real

    res.set(0, 3*nstates+0, ave_en_se)  # Average SE energy
    res.set(0, 3*nstates+1, ave_en_sh)  # Average SH energy
    res.set(0, 3*nstates+2, tot_se)     # Total SE population
    res.set(0, 3*nstates+3, tot_sh)     # Total SH population

    #print "End of traj_statistics...\n"

    return res




def printout(t, res, outfile):
    """
    t        [float]       - time [a.u.] 
    res      [MATRIX(1,N)] - information to be printed out
    outfile  [string]      - filename where we'll print everything out

    """

    N = res.num_of_cols

    line = "%8.5f " % (t)
    for i in xrange(N):
        line = line + " %8.5f " % (res.get(0,i))
    line = line + " \n"

    f = open(outfile, "a") 
    f.write(line)
    f.close()




def run(H_vib, params):
    """
    The main procedure to run NA-MD calculations according to the PYXAID2 workflow

    H_vib                     [list of lists of CMATRIX] - the vibronic Hamiltonian for all data sets and all time-points

    H_vib[idata][istep].get(i,j) - i,j matrix element for the data set ```idata``` and step in that data set ```istep```
    

    === Required parameter keys: ===

    params["T"]                [double] - temperature of nuclear/electronic dynamics [in K, default: 300.0]
    params["ntraj"]            [int] - the number of stochastic surface hopping trajectories [default: 1]
    params["sh_method"]        [int] - how to compute surface hopping probabilities [default: 1]
                               Options:
                               0 - MSSH
                               1 - FSSH
    params["decoherence_method"]  [int] - selection of decoherence method [default: 0]
                               Options:
                               0 - no decoherence
                               1 - ID-A
                               2 - MSDM
                               3 - DISH
    params["dt"]               [double] - nuclear dynamics integration time step [in a.u. of time, default: 41.0]
    params["nsteps"]           [int] - the length of the NA-MD trajectory
    params["Boltz_opt"]        [int] - option to select a probability of hopping acceptance [default: 3]
                               Options:
                               0 - all proposed hops are accepted - no rejection based on energies
                               1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
                               2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
                               3 - generalization of "1", but actually it should be changed in case there are many degenerate levels
    params["istate"]           [int] - index of the initial spin-adapted state [default: 0]
    params["init_times"]       [list of ints] - indices of the starting point in the provided data arrays [default: [0]]
    params["outfile"]          [string] - the name of the file where to print populations and energies of states [default: "_out.txt"]    


    === Required by the get_Hvib2() ===

    params["nstates"]          [int] - how many lines/columns in the file - the total number of states set in step3
    params["nfiles"]           [int] - hom many files there are to read for each data set
    params["data_set_paths"]   [list of strings] - define the pathes of the directories where the vibronic Hamiltonians for
                               different data sets (independent MD trajectories) are located
    params["Hvib_re_prefix"]   [string] - prefixes of the files with real part of the vibronic Hamiltonian at time t
    params["Hvib_re_suffix"]   [string] - suffixes of the files with real part of the vibronic Hamiltonian at time t
    params["Hvib_im_prefix"]   [string] - prefixes of the files with imaginary part of the vibronic Hamiltonian at time t
    params["Hvib_im_suffix"]   [string] - suffixes of the files with imaginary part of the vibronic Hamiltonian at time t

   
    ```The full name of the vibronic Hamiltonian files will be:
    
    params["data_set_path"]+params["Hvib_re_prefix"]+integer(time step)+params["Hvib_re_suffix"] - for real part
    params["data_set_path"]+params["Hvib_im_prefix"]+integer(time step)+params["Hvib_im_suffix"] - for imaginary part

    Example: 
    Say, the directory "/home/alexeyak/test/step3/res0" contains files:
    Hvib_0_re, Hvib_1_re, .... ,    Hvib_999_re
    Hvib_0_im, Hvib_1_im, .... ,    Hvib_999_im

    Then set:
    params["data_set_paths"] = ["/home/alexeyak/test/step3/res0/"]
    params["Hvib_re_prefix"] = "Hvib_"
    params["Hvib_re_suffix"] = "_re"
    params["Hvib_im_prefix"] = "Hvib_"
    params["Hvib_im_suffix"] = "_im"
    ```
     
    """

    critical_params = [ "nstates", "nsteps" ] 
    default_params = { "T":300.0, "ntraj":1, 
                       "sh_method":1, "decoherence_method":0, "dt":41.0, "Boltz_opt":3,
                       "istate":0, "init_times":[0], "outfile":"_out.txt",
                       "on-the-fly-qsh":False }
    comn.check_input(params, default_params, critical_params)


    rnd = Random()

#    if params["on-the-fly-qsh"]:
        # To ensure consistency, we'll borrow some of the parameters from the
        # QSH calculations
#        params.update(params["qsh-params"])
        


    nsteps = params["nsteps"]
    nstates = params["nstates"]
    ndata = len(params["data_set_paths"])
    ntraj = params["ntraj"]
    nitimes = len(params["init_times"])
    Ntraj = ndata * nitimes * ntraj

    T = params["T"]
    bolt_opt = params["Boltz_opt"]
    dt = params["dt"]

    res = MATRIX(nsteps, 3*nstates+5)


    #========== Compute PARAMETERS  ===============
    # Decoherence times for DISH
    tau, decoh_rates = dectim.decoherence_times_ave(H_vib, params["init_times"], nsteps, 1) 

    #========== Initialize the DYNAMICAL VARIABLES  ===============
    # TD-SE coefficients and active state indices
    Coeff, istate = [], []

    # Coherence times and coherence intervals for DISH
    t_m, tau_m = [], []

    for tr in xrange(Ntraj):
        istate.append(params["istate"])
        Coeff.append(CMATRIX(nstates, 1)); 
        Coeff[tr].set(params["istate"], 1.0, 0.0)
        t_m.append(MATRIX(nstates,1))
        tau_m.append(MATRIX(nstates,1))  

           
    # Prepare the output file
    f = open(params["outfile"],"w"); f.close()

    #=============== Entering the DYNAMICS ========================
    for i in xrange(nsteps):  # over all evolution times


        #============== Analysis of the Dynamics  =================
        # Compute the averages
        res_i = traj_statistics(i, Coeff, istate, H_vib, params["init_times"])

        # Print out into a file
        printout(i*dt, res_i, params["outfile"])

        # Update the overal results matrix
        res.set(i,0, i*dt)
        push_submatrix(res, res_i, Py2Cpp_int([i]), Py2Cpp_int(range(1,3*nstates+5)) )


        #=============== Propagation ==============================
        for idata in xrange(ndata):   # over all MD trajectories (data sets)

            for it_indx in xrange(nitimes): # over all initial times

                it = params["init_times"][it_indx]

                for tr in xrange(ntraj):  # over all stochastic trajectories

                    Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                    #============== Propagation: TD-SE and surface hopping ==========
        
                    # Coherent evolution amplitudes
                    propagate_electronic(dt, Coeff[Tr], H_vib[idata][it+i])   # propagate the electronic DOFs

        
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
                            Coeff[Tr] = msdm(Coeff[Tr], dt, istate[Tr], decoh_rates)
                    
                        istate[Tr], Coeff[Tr] = tsh.hopping(Coeff[Tr], H_vib[idata][it+i], istate[Tr], params["sh_method"], do_collapse, ksi, ksi2, dt, T, bolt_opt)
                    
                    elif params["decoherence_method"] in [3]:  # DISH
                    
                        tau_m[Tr] = coherence_intervals(Coeff[Tr], decoh_rates)
                        istate[Tr] = tsh.dish_py(Coeff[Tr], istate[Tr], t_m[Tr], tau_m[Tr], H_vib[idata][it+i], bolt_opt, T, ksi, ksi2)
                        t_m[Tr] += dt

        
    return res


