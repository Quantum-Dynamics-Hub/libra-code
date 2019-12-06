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

  Implementation of the Landau-Zener transition probability.

"""

import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

#import libra_py.common_utils as comn
import util.libutil as comn

import libra_py.units as units
import libra_py.probabilities as prob
import libra_py.tsh as tsh
from . import decoherence_times
from . import step4



def Belyaev_Lebedev(Hvib, params):
    """
    Computes the Landau-Zener hopping probabilities based on the energy levels
    according to: 
    (1) Belyaev, A. K.; Lebedev, O. V. Phys. Rev. A, 2011, 84, 014701

    See also:
    (2) Xie, W.; Domcke, W. J. Chem. Phys. 2017, 147, 184114
    (3) Crespo-Otero, R.; Barbatti, M. Chem. Rev. 2018, 188, 7026 - section: 3.2.3

    Specifics:
    1) The estimation of d^2E_ij / dt^2 is based on the 3-point Lagrange interpolation
    2) This is done within the NBRA

    Args:
        Hvib (list of CMATRIX(nstates,nstates) ):  vibronic Hamiltonians along the trajectory
        params ( dictionary ): control parameters

            * **params["dt"]** ( double ): time distance between the adjacent data points [ units: a.u., defaut: 41.0 ]
            * **params["T"]** ( double ): temperature of the nuclear sub-system [ units: K, default: 300.0 ]
            * **params["Boltz_opt_BL"]** ( int ): option to select a probability of hopping acceptance [default: 1]                
                Options:

                    - 0 - all proposed hops are accepted - no rejection based on energies
                    - 1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
                    - 2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
                    - 3 - generalization of "1", but actually it should be changed in case there are many degenerate levels

            * **params["gap_min_exception"]** ( int ): option to handle the situation when extrapolated gap minimum is negative
                Options:

                    - 0 - set to zero [ default ]
                    - 1 - use the mid-point gap

            * **params["target_space"]** ( int ): how to select the space of target states for each source state
                Options:

                    - 0 - only adjacent states 
                    - 1 - all states available [ default ]

     
    """

    # Control parameters
    critical_params = [  ]
    default_params = { "T":300.0, "Boltz_opt_BL":1, "dt":41.0, "gap_min_exception":0, "target_space":1 }
    comn.check_input(params, default_params, critical_params)

    boltz_opt = params["Boltz_opt_BL"]
    T = params["T"]
    dt = params["dt"]
    gap_min_exception = params["gap_min_exception"]
    target_space = params["target_space"]


    # Data dimensions 
    nsteps = len(Hvib)
    nstates= Hvib[0].num_of_cols


    # Pre-compute the energy gaps along the trajectory 
    dE = decoherence_times.energy_gaps(Hvib)


    """
    Compute the probabilities based on the LZ formula
    P(i,j) - the convention is: the probability to go from j to i
    This will make the Markov state propagation more convenient
    """


    P = []
    P.append( MATRIX(nstates, nstates) )
    for i in range(0,nstates):    
        P[0].set(i,i, 1.0)

    for n in range(1, nsteps-1):
        P.append(MATRIX(nstates, nstates))
 

        # Belyaev-Lebedev probabilities
        # Find the minima of the |E_i - E_j| for all pair of i and j    
        for j in range(0, nstates):      # source 


            # By doing this, we'll consider the transitions to only adjacent states
            targets = []

            if target_space == 0:
                # Target states are only the states adjacent to the source state

                if j == 0:
                    targets = [1]
                elif j == nstates - 1:
                    targets = [nstates - 2]
                else:
                    targets = [j-1, j+1]

            elif target_space == 1:
                # All states can be the potential targets

                targets = range(0, nstates)




            # with how many other states, does the current (source) state has minima?
            # apparently, this can not be more than 2
            normalization = 0.0
            
            #for i in targets:   # targets 
            for i in targets:

                # Interpolation is based on the 3-points Lagrange interpolant
                # http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html 

                if (dE[n-1].get(i,j)>dE[n].get(i,j) and dE[n].get(i,j)<dE[n+1].get(i,j)):

                                        
                    denom = dE[n-1].get(i,j) - 2.0*dE[n].get(i,j) + dE[n+1].get(i,j)                  
                    if denom > 0.0:
                       
                        t_min = 0.5*(dE[n-1].get(i,j) - dE[n+1].get(i,j))*dt/denom

                        if t_min<-dt or t_min > dt:
                            print("Error determining t_min in the interpolation!\n")
                            print("Exiting...\n")
                            sys.exit(0)

                        gap_min = 0.5*(t_min*(t_min - dt)*dE[n-1].get(i,j) 
                          - 2.0*(t_min + dt)*(t_min - dt)*dE[n].get(i,j) + 
                                       t_min*(t_min + dt)*dE[n+1].get(i,j) )/(dt*dt)

                        if gap_min<0.0:
                            if gap_min_exception==0:
                                gap_min = 0.0
                            elif gap_min_exception==1:
                                gap_min = dE[n].get(i,j)

                        if gap_min > dE[n-1].get(i,j) or gap_min > dE[n+1].get(i,j):
                            print("Error: the extrapolated gap is larger than the bounding values!\n")
                            print("Exiting...\n")
                            sys.exit(0)

                        second_deriv = denom/(dt*dt)

                        argg = (gap_min**3) / second_deriv

                        p = math.exp(-0.5*math.pi*math.sqrt(argg) )


                        # Optionally, can correct transition probabilitieis to 
                        # account for Boltzmann factor
                        if i!=j:

                            E_new = Hvib[n].get(i,i).real  # target
                            E_old = Hvib[n].get(j,j).real  # source

                            bf = 1.0
                            if E_new > E_old:
                                # Notice how we use gap_min rather than E_new - E_old in this case
                                bf = tsh.boltz_factor(gap_min, 0.0, T, boltz_opt)

                            if bf>1.0:
                                print("Error: Boltzmann scaling factor can not be larger 1.0 = ",bf)
                                sys.exit(0)


                            P[n].set(i,j, p*bf)      # Probability to go j->i
                            normalization = normalization + p*bf
                                 

            if normalization<1.0:
                P[n].set(j,j, 1.0 - normalization)
            else:
                P[n].add(j,j, 0.0)

                scl = 1.0/normalization
                P[n].scale(-1, j, scl)
 
               
       

    P.append( MATRIX(nstates, nstates) )
    for i in range(0,nstates):    
        P[nsteps-1].set(i,i, 1.0)

            
    return P



def run(H_vib, params):
    """
    Main function to run the SH calculations based on the Landau-Zener hopping
    probabilities, all within the NBRA. The probabilities are implemented according to
    the Belyaev-Lebedev work.

    Args:
        params ( dictionary ): control parameters

        * **params["dt"]** ( double ) : time distance between the adjacent data points [ units: a.u.; default: 41.0 a.u.]
        * **params["ntraj"]** ( int ) : how many stochastic trajectories to use in the ensemble [ defult: 1]
        * **params["nsteps"]** ( int ) : how nuclear steps in the trajectory to be computed [ defult: 1]
        * **params["istate"]** ( int ) : index of the starting state (within those used in the active_space) [ default: 0]
        * **params["Boltz_opt"]** ( int ) : option to control the acceptance of the proposed hops 
            Options:

            - 0 - all proposed hops are accepted - no rejection based on energies
            - 1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach) [ default ]
            - 2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
            - 3 - generalization of "1", but actually it should be changed in case there are many degenerate levels

        * **params["Boltz_opt_BL"]** ( int ) : what type of hop acceptance scheme to incorporate into the BL probabilities 
            Options:

            - 0 - all proposed hops are accepted - no rejection based on energies [ default ]
            - 1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach) 
            - 2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
            - 3 - generalization of "1", but actually it should be changed in case there are many degenerate levels
        
        * **params["T"]** ( double ) : temperature of the nuclei - affects the acceptance probabilities [ units: K, default: 300.0 K]
        * **params["do_output"]** ( Boolean ) : wheather to print out the results into a file [ default: True ]
        * **params["outfile"]** ( string ) : the name of the file, where all the results will be printed out [ default: "_out.txt" ]
        * **params["do_return"]** ( Boolean ) : wheather to construct the big matrix with all the result [ default: True ]
        * **params["evolve_Markov"]** ( Boolean ) : wheather to propagate the "SE" populations via Markov chain [ default: True ]
        * **params["evolve_TSH"]** ( Boolean ) : wheather to propagate the "SH" populations via TSH with many trajectories [ default: True ]

    """


    critical_params = [  ]
    default_params = { "dt":41.0, "ntraj":1, "nsteps":1, "istate":0, 
                       "Boltz_opt":1, "Boltz_opt_BL":1, "T":300.0,
                       "do_output":True, "outfile":"_out.txt", "do_return":True,
                       "evolve_Markov":True, "evolve_TSH":True }
    comn.check_input(params, default_params, critical_params)
    
    rnd = Random()

    ndata = len(H_vib)
    nsteps = params["nsteps"]
    nstates= H_vib[0][0].num_of_cols
    dt = params["dt"]
    do_output = params["do_output"]
    do_return = params["do_return"]
    ntraj = params["ntraj"]
    boltz_opt = params["Boltz_opt"]
    T = params["T"]
    evolve_Markov = params["evolve_Markov"]
    evolve_TSH = params["evolve_TSH"]


    res = MATRIX(nsteps, 3*nstates+5)

    #===== Precompute hopping probabilities ===
    P = []
    itimes = params["init_times"]
    nitimes = len(itimes)

    for idata in range(0,ndata):
        p = Belyaev_Lebedev(H_vib[idata], params)
        P.append(p)


    #========== Initialize the DYNAMICAL VARIABLES  ===============
    # State populations and active state indices
    Pop, istate = [], []
    Ntraj = ndata * nitimes * ntraj

    for tr in range(0,Ntraj):
        istate.append(params["istate"])
        Pop.append(CMATRIX(nstates, 1)); 
        Pop[tr].set(params["istate"], 1.0, 0.0)


    #=============== Entering the DYNAMICS ========================
    for i in range(0,nsteps):  # over all evolution times

        #============== Analysis of the Dynamics  =================
        # Compute the averages
        #res_i = step4.traj_statistics(i, Coeff, istate, H_vib, itimes)
        #res_i = step4.traj_statistics2(i, Pop, istate, H_vib, itimes)
        res_i = step4.traj_statistics2_fast(i, Pop, istate, H_vib, itimes)

        # Print out into a file
        if do_output==True:
            step4.printout(i*dt, res_i, params["outfile"])

        # Update the overal results matrix
        res.set(i,0, i*dt)
        if do_return==True:
            push_submatrix(res, res_i, Py2Cpp_int(list([i])), Py2Cpp_int(list(range(1,3*nstates+5))) )

        #=============== Propagation ==============================
        for idata in range(0,ndata):   # over all data sets (MD trajectories)

            for it_indx in range(0,nitimes): # over all initial times within each MD trajectory

                it = itimes[it_indx]

                for tr in range(0,ntraj):  # over all stochastic trajectories

                    Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                    #============== Propagation: TD-SE and surface hopping ==========
        
                    # Evolve the Markov process.
                    # The convention is:
                    # P(i,j) - the probability to go from j to i
                    if evolve_Markov==True:
                        Pop[Tr] = CMATRIX(P[idata][i]) * Pop[Tr]


                    if evolve_TSH==True:        

                        # Surface hopping 
                        ksi  = rnd.uniform(0.0, 1.0)
                        
                        # Proposed hop:
                        st_new = tsh.hop_py(istate[Tr], P[idata][i].T(), ksi)  
                        
                        # Accept the proposed hop with the Boltzmann probability
                        E_new = H_vib[idata][i].get(st_new,st_new).real
                        E_old = H_vib[idata][i].get(istate[Tr], istate[Tr]).real
                        de = E_new - E_old
                        
                        if de>0.0:
                            bf = tsh.boltz_factor(E_new, E_old, T, boltz_opt)
                            ksi  = rnd.uniform(0.0, 1.0)
                            if ksi < bf:
                                istate[Tr] = st_new                  
                        else:
                            istate[Tr] = st_new
        
    return res


