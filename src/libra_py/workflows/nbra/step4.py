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
#
#  
#
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


def get_Hvib(params):
    """
    Read a single set of generic vibronic Hamiltonian files 

    Required parameter keys:

    params["nstates"]          [int] - how many lines/columns in the file 
    params["nsteps"]           [int] - how many files to read, starting from index 0
    params["Hvib_re_prefix"]   [string] - prefixes of the files with real part of the MO overlaps at time t
    params["Hvib_re_suffix"]   [string] - suffixes of the files with real part of the MO overlaps at time t
    params["Hvib_im_prefix"]   [string] - prefixes of the files with imaginary part of the MO overlaps at time t
    params["Hvib_im_suffix"]   [string] - suffixes of the files with imaginary part of the MO overlaps at time t

    """

    critical_params = ["nstates", "nsteps", "Hvib_re_prefix", "Hvib_im_prefix"]
    default_params = { "Hvib_re_suffix":"_re", "Hvib_im_suffix":"_im"}
    comn.check_input(params, default_params, critical_params)

    nsteps = params["nsteps"]
    nstates = params["nstates"]  # the number of states in the input files
    active_space = range(nstates)

    Hvib = []

    for i in range(0,nsteps):

        filename_re = params["Hvib_re_prefix"]+str(i)+params["Hvib_re_suffix"]
        filename_im = params["Hvib_im_prefix"]+str(i)+params["Hvib_im_suffix"]
        hvib = comn.get_matrix(nstates, nstates, filename_re, filename_im, active_space )
        Hvib.append(hvib)

    return Hvib



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

    Ntraj = len(Coeff)               # total number of trajectpry = data set size x number of init times x number of SH trajectories
    nstates = Coeff[0].num_of_rows   # the number of states
    ndata = len(Hvib)                # how many data sets
    nitimes = len(itimes)            # how many initial times
    ntraj = Ntraj/(ndata * nitimes)  # how many stochastic SH trajectories per data set/initial condition


    # Update SE-derived density matrices
    denmat_se = tsh.amplitudes2denmat(Coeff)   # list of Ntraj CMATRIX elements 

    # Use SE-derived density matrices to make SH-derived density matrices
    denmat_sh = []
    H_vib = []
    H_vib_ave = CMATRIX(nstates,nstates)  # Hvib averaged over the data sets/initial times

    for idata in xrange(ndata):
        for it_indx in xrange(nitimes):
            it = itimes[it_indx]

            H_vib_ave += Hvib[idata][it+i]

            for tr in xrange(ntraj):

                Tr = idata*ndata + it_indx*nitimes + tr

                denmat_sh.append(CMATRIX(nstates, nstates))
                denmat_sh[Tr].set(istate[Tr],istate[Tr], 1.0, 0.0)
                H_vib.append(Hvib[idata][it+i])   

    H_vib_ave *= (1.0/float(ndata * nitimes))

    # Update TSH-ensemble-averaged SE and SH populations 
    ave_pop_sh, ave_pop_se = tsh.ave_pop(denmat_sh, denmat_se)
    ave_en_sh,  ave_en_se  = tsh.ave_en(denmat_sh, denmat_se, H_vib)

    # Save the computed data into a matrix to be output
    res = MATRIX(1, 3*nstates+4) 
   
    tot_sh, tot_se = 0.0, 0.0
    for j in xrange(nstates):
        res.set(i, 3*j+0, H_vib_ave.get(j,j).real)   # Energy of the state j
        res.set(i, 3*j+1, ave_pop_se.get(j,j).real)  # SE population
        res.set(i, 3*j+2, ave_pop_sh.get(j,j).real)  # SH population

        tot_se += ave_pop_se.get(j,j).real
        tot_sh += ave_pop_sh.get(j,j).real

    res.set(i, 3*nstates+0, ave_en_se)  # Average SE energy
    res.set(i, 3*nstates+1, ave_en_sh)  # Average SH energy
    res.set(i, 3*nstates+2, tot_se)     # Total SE population
    res.set(i, 3*nstates+3, tot_sh)     # Total SH population

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




def run(params):
    """
    The main procedure to run NA-MD calculations according to the PYXAID2 workflow

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

    params["Boltz_opt"]        [int] - option to select a probability of hopping acceptance [default: 3]
                               Options:
                               0 - all proposed hops are accepted - no rejection based on energies
                               1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
                               2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
                               3 - generalization of "1", but actually it should be changed in case there are many degenerate levels
    params["istate"]           [int] - index of the initial spin-adapted state [default: 0]
    params["init_times"]       [list of ints] - indices of the starting point in the provided data arrays [default: [0]]
    params["outfile"]          [string] - the name of the file where to print populations and energies of states [default: "_out.txt"]    


    === Required by the get_data() ===

    params["nstates"]          [int] - how many lines/columns in the file - the total number of spin-orbitals
    params["nsteps"]           [int] - how many files to read, starting from index 0
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

    critical_params = [ "nstates", "nsteps", "data_set_paths" ]
    default_params = { "T":300.0, "ntraj":1, 
                       "sh_method":1, "decoherence_method":0, "dt":41.0, "Boltz_opt":3,
                       "istate":0, "init_times":[0], "outfile":"_out.txt" }
    comn.check_input(params, default_params, critical_params)


    rnd = Random()


    nsteps = params["nsteps"]
    nstates = params["nstates"]
    ntraj = params["ntraj"]
    ndata = len(params["data_set_paths"])
    nitimes = len(params["init_times"])
    Ntraj = ndata * nitimes * ntraj

    T = params["T"]
    bolt_opt = params["Boltz_opt"]
    dt = params["dt"]
    do_collapse = 0
    if params["decoherence_method"]==1:
        do_collapse = 1


    res = MATRIX(nsteps, 3*nstates+5)

    #======== Read in the vibronic Hamiltonian along the trajectory for each data set ============
    H_vib = []

    for idata in params["data_set_paths"]:   # over all MD trajectories (data sets)
        prms = dict(params)    
        prms.update({"Hvib_re_prefix": idata+params["Hvib_re_prefix"] })
        prms.update({"Hvib_im_prefix": idata+params["Hvib_im_prefix"] })                

        h_vib = get_Hvib(prms)  
        H_vib.append(h_vib)


    #========== Compute PARAMETERS  ===============
    # Decoherence times for DISH
    tau, decoh_rates = dectim.decoherence_times(H_vib, params["init_times"], nsteps, 1) 

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

           
    #=============== Entering the DYNAMICS ========================
    for i in xrange(nsteps):  # over all evolution times


        #============== Analysis of the Dynamics  =================
        # Compute the averages
        res_i = traj_statistics(i, Coeff, istate, H_vib, params["init_times"])

        # Print out into a file
        printout(i*dt, res_i, params["outfile"])

        # Update the overal results matrix
        res.set(i,0, i*dt)
        push_submatrix(res, res_i, Py2Cpp_int([i]), Py2Cpp_int(range(1,3*nstates+6)) )


        #=============== Propagation ==============================
        for idata in xrange(ndata):   # over all MD trajectories (data sets)

            for it_indx in xrange(nitimes): # over all initial times

                it = params["init_times"][it_indx]

                for tr in xrange(ntraj):  # over all stochastic trajectories

                    Tr = idata*ndata + it_indx*nitimes + tr

                    #============== Propagation: TD-SE and surface hopping ==========
        
                    # Coherent evolution amplitudes
                    propagate_electronic(dt, Coeff[Tr], H_vib[idata][it+i]])   # propagate the electronic DOFs
        
                    # Surface hopping 
                    ksi  = rnd.uniform(0.0, 1.0)
                    ksi2 = rnd.uniform(0.0, 1.0)
        
                    if params["decoherence_method"] in [0, 1, 2]:
                    
                        if params["decoherence_method"]==0:    # No decoherence
                            pass
                        elif params["decoherence_method"]==1:  # ID-A, taken care of in the tsh.hopping
                            pass
                        elif params["decoherence_method"]==2:  # MSDM
                            Coeff[Tr] = msdm(Coeff[Tr], dt, istate[Tr], decoh_rates)
                    
                        istate[Tr] = tsh.hopping(Coeff[Tr], H_vib[idata][it+i], istate[Tr], params["sh_method"], do_collapse, ksi, ksi2, dt, T, bolt_opt)
                    
                    elif params["decoherence_method"] in [3]:  # DISH
                    
                        tau_m[Tr] = coherence_intervals(Coeff[Tr], decoh_rates)
                        istate[Tr] = tsh.dish_py(Coeff[Tr], istate[Tr], t_m[Tr], tau_m[Tr], H_vib[idata][it+i], bolt_opt, T, ksi, ksi2)
                        t_m[Tr] += dt
        
    return res


