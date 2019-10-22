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
            * **params["Boltz_opt"]** ( int ): option to select a probability of hopping acceptance [default: 3]
                Options:

                - 0 - all proposed hops are accepted - no rejection based on energies
                - 1 - proposed hops are accepted with exp(-E/kT) probability - the old (hence the default approach)
                - 2 - proposed hops are accepted with the probability derived from Maxwell-Boltzmann distribution - more rigorous
                - 3 - generalization of "1", but actually it should be changed in case there are many degenerate levels

     
    """

    # Control parameters
    critical_params = [  ]
    default_params = { "T":300.0, "Boltz_opt":3, "dt":41.0 }
    comn.check_input(params, default_params, critical_params)

    boltz_opt = params["Boltz_opt"]
    T = params["T"]
    dt = params["dt"]


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
        for i in range(0,nstates):             # target
            for j in range(i+1, nstates):      # source 

                # Interpolation is based on the 3-points Lagrange interpolant
                # http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html 
                p = 0.0
                if (dE[n-1].get(i,j)>dE[n].get(i,j) and dE[n].get(i,j)<dE[n+1].get(i,j)):
                    
                    denom = dE[n-1].get(i,j) - 2.0*dE[n].get(i,j) + dE[n+1].get(i,j) 
                    if denom > 0.0:
                        argg = (dE[n].get(i,j)**3) / denom
                        p = math.exp(-0.5*math.pi*dt*math.sqrt(argg) )
                else:
                    p = 0.0   # no transitions is not a minimum

                P[n].set(i,j, p)
                P[n].set(j,i, p)

        
        # Optionally, can correct transition probabilitieis to 
        # account for Boltzmann factor
        for i in range(0,nstates):        # target
            for j in range(0,nstates):    # source 

                if i!=j:

                    E_new = Hvib[n].get(i,i).real  # target
                    E_old = Hvib[n].get(j,j).real  # source
                    bf = 1.0
                    if E_new > E_old:
                        bf = tsh.boltz_factor(E_new, E_old, T, boltz_opt)
                        if bf>1.0:
                            print("Error: Boltzmann scaling factor can not be larger 1.0 = ",bf)
                            #sys.exit(0)
                        P[n].scale(i,j, bf)


        # Compute the probability of staying on the same state j (source)
        # P(j,j) = 1 - sum_(i!=j) { P(i,j) }
        #
        # The convention is:
        # P(i,j) - the probability to go from j to i

        for j in range(0,nstates): # for all source states            

            tot = 0.0  # Total probability to leave state j
            for i in range(0,nstates):     # all target states
                if i!=j:                  # but j
                    tot += P[n].get(i,j)

            # Compute the probability to stay on state j
            P[n].set(j,j, 1.0 - tot)


    P.append( MATRIX(nstates, nstates) )
    for i in range(0,nstates):    
        P[nsteps-1].set(i,i, 1.0)



            
    return P


#        # Normalize probabilities if needed
#        if tot>1.0:
#            for i in xrange(nstates):        
#            P[n].scale(i,j, (1.0/tot))
#            tot = 1.0


def run(H_vib, params):
    """
    Main function to run the SH calculations based on the Landau-Zener hopping
    probabilities, all within the NBRA. The probabilities are implemented according to
    the Belyaev-Lebedev work.

    ===== Modeling params ===== 

    params["dt"]             [double, a.u.] - time distance between the adjacent data points
    params["ntraj"]          [int] - how many stochastic trajectories to use in the ensemble
    params["istate"]         [int] - index of the starting state (within those used in the active_space - see above)
    params["do_output"]      [string] - wheather to print out the results into a file
    params["outfile"]        [string] - the name of the file, where all the results will be printed out
    params["T"]              [double, K] - temperature of the simulation

    """


    critical_params = [  ]
    default_params = { "T":300.0, "ntraj":1000, "nsteps":1,"istate":0, 
                       "sh_method":1, "decoherence_method":0, "dt":41.0, 
                       "Boltz_opt":3,
                       "do_output":False, "outfile":"_out.txt" }
    comn.check_input(params, default_params, critical_params)
    
    rnd = Random()

    ndata = len(H_vib)
    nsteps = params["nsteps"]
    nstates= H_vib[0][0].num_of_cols
    dt = params["dt"]
    do_output = params["do_output"]
    ntraj = params["ntraj"]
    boltz_opt = params["Boltz_opt"]
    T = params["T"]


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
        res_i = step4.traj_statistics2(i, Pop, istate, H_vib, itimes)

        # Print out into a file
        step4.printout(i*dt, res_i, params["outfile"])

        # Update the overal results matrix
        res.set(i,0, i*dt)
        push_submatrix(res, res_i, Py2Cpp_int(list([i])), Py2Cpp_int(list(range(1,3*nstates+5))) )


        #=============== Propagation ==============================
        for idata in range(0,ndata):   # over all data sets (MD trajectories)

            for it_indx in range(0,nitimes): # over all initial times within each MD trajectory

                it = itimes[it_indx]

                for tr in range(0,ntraj):  # over all stochastic trajectories

                    Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                    #============== Propagation: TD-SE and surface hopping ==========
        
                    # Evolve Markov process.
                    # The convention is:
                    # P(i,j) - the probability to go from j to i
                    Pop[Tr] = CMATRIX(P[idata][i]) * Pop[Tr]

        
                    # Surface hopping 
                    ksi  = rnd.uniform(0.0, 1.0)
                    ksi1 = rnd.uniform(0.0, 1.0)

                    
                    # Proposed hop:
                    st_new = tsh.hop_py(istate[Tr], P[idata][i].T(), ksi)  

                    # Accept the proposed hop with the Boltzmann probability
                    E_new = H_vib[idata][i].get(st_new,st_new).real
                    E_old = H_vib[idata][i].get(istate[Tr], istate[Tr]).real
                    de = E_new - E_old
                    
                    if de>0.0:
                        bf = tsh.boltz_factor(E_new, E_old, T, boltz_opt)
                        if ksi1 < bf:
                            istate[Tr] = st_new                  
                    else:
                        istate[Tr] = st_new
        
    return res


