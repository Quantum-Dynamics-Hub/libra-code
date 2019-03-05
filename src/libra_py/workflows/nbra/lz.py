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

import libra_py.common_utils as comn
import libra_py.units as units
import libra_py.probabilities as prob
import libra_py.tsh as tsh
import decoherence_times
import step4



def Belyaev_Lebedev(Hvib, dt):
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

    Parameters:
    Hvib [list of CMATRIX] - vibronic Hamiltonians along the trajectory
    dt [double, a.u.] - time distance between the adjacent data points
     
    """

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
    for i in xrange(nstates):    
        P[0].set(i,i, 1.0)

    for n in xrange(1, nsteps-1):
        # Find the minima of the |E_i - E_j| for all pair of i and j    
        P.append(MATRIX(nstates, nstates))
        for i in xrange(nstates):             # target
            for j in xrange(i+1, nstates):    # source 

                # Interpolation is based on the 3-points Lagrange interpolant
                # http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html 
                p = 0.0
                if (dE[n-1].get(i,j)>dE[n].get(i,j) and dE[n].get(i,j)<dE[n+1].get(i,j)):
                    
                    denom = dE[n-1].get(i,j) - 2.0*dE[n].get(i,j) + dE[n+1].get(i,j) 
                    if denom > 0.0:
                        argg = (dE[n].get(i,j)**3) / denom
                        p = math.exp(-2.0*math.pi*dt*math.sqrt(argg) )
                else:
                    p = 0.0

                P[n].set(i,j, p)
                P[n].set(j,i, p)


        # Compute the probability of staying on the same state
        for j in xrange(nstates):         # source
            # Compute total probabilities to leave state j
            tot = 0.0
            for i in xrange(nstates):     # target
                if i!=j:
                    tot += P[n].get(i,j)

            # Normalize probabilities if needed
            if tot>1.0:
                for i in xrange(nstates):        
                    P[n].scale(i,j, (1.0/tot))
                tot = 1.0

            # Compute the probability to stay on state j
            P[n].set(j,j, 1.0 - tot)


    P.append( MATRIX(nstates, nstates) )
    for i in xrange(nstates):    
        P[nsteps-1].set(i,i, 1.0)

            
    return P


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
    nsteps = len(H_vib[0])
    nstates= H_vib[0][0].num_of_cols
    dt = params["dt"]
    do_output = params["do_output"]
    ntraj = params["ntraj"]
    boltz_opt = params["Boltz_opt"]
    T = params["T"]


    res = MATRIX(nsteps, 3*nstates+5)

    #===== Precompute hopping probabilities ===
    P = []
    itimes = [0]
    nitimes = len(itimes)

    for idata in xrange(ndata):
        p = Belyaev_Lebedev(H_vib[idata], dt)
        P.append(p)


    #========== Initialize the DYNAMICAL VARIABLES  ===============
    # TD-SE coefficients and active state indices
    Coeff, istate = [], []

    for tr in xrange(ntraj):
        istate.append(params["istate"])
        Coeff.append(CMATRIX(nstates, 1)); 
        Coeff[tr].set(params["istate"], 1.0, 0.0)


    #=============== Entering the DYNAMICS ========================
    for i in xrange(nsteps):  # over all evolution times

        #============== Analysis of the Dynamics  =================
        # Compute the averages
        res_i = step4.traj_statistics(i, Coeff, istate, H_vib, itimes)

        # Print out into a file
        step4.printout(i*dt, res_i, params["outfile"])

        # Update the overal results matrix
        res.set(i,0, i*dt)
        push_submatrix(res, res_i, Py2Cpp_int([i]), Py2Cpp_int(range(1,3*nstates+5)) )


        #=============== Propagation ==============================
        for idata in xrange(ndata):   # over all MD trajectories (data sets)

            for it_indx in xrange(nitimes): # over all initial times

                it = itimes[it_indx]

                for tr in xrange(ntraj):  # over all stochastic trajectories

                    Tr = idata*(nitimes*ntraj) + it_indx*(ntraj) + tr

                    #============== Propagation: TD-SE and surface hopping ==========
        
                    # Coherent evolution amplitudes
                    #propagate_electronic(dt, Coeff[Tr], H_vib[idata][it+i])   # propagate the electronic DOFs
                    Coeff[Tr] = CMATRIX(P[idata][i].T()) * Coeff[Tr]

        
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


