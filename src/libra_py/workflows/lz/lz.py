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

import libra_py.workflows.common_utils as comn
import libra_py.units as units
import libra_py.probabilities as prob
import libra_py.tsh as tsh



def get_data(params):
    """
    Read in all the Hvib files (we actually need only energies)

    Required parameter keys:

    params["norbitals"]      [int] - how many lines/columns in the file
    params["active_space"]   [list of ints] - which orbitals we care about (indexing starts with 0)
    params["Hvib_re_prefix"] [string] - prefixes of the files with real part of the Hamiltonian
    params["Hvib_im_prefix"] [string] - prefixes of the files with imaginary part of the Hamiltonian
    params["Hvib_re_suffix"] [string] - suffixes of the files with real part of the Hamiltonian
    params["Hvib_im_suffix"] [string] - suffixes of the files with imaginary part of the Hamiltonian
    params["nsteps"]         [int] - how many files to read

    """

    critical_params = ["norbitals", "active_space", "Hvib_re_prefix", "Hvib_im_prefix", "nsteps" ]
    default_params = { "Hvib_re_suffix":"_re", "Hvib_im_suffix":"_im"}
    comn.check_input(params, default_params, critical_params)


    norbitals = params["norbitals"]  # the number of orbitals in the input files
    active_space = params["active_space"]
    nstates = len(active_space)
    nsteps = params["nsteps"]  # how many steps 
    
    Hvib = []    
    for step in xrange(0, nsteps): # how many files we have
        filename_re = params["Hvib_re_prefix"]+str(step)+params["Hvib_re_suffix"]
        filename_im = params["Hvib_im_prefix"]+str(step)+params["Hvib_im_suffix"]
        hvib = comn.get_matrix(norbitals, norbitals, filename_re, filename_im, active_space)
        Hvib.append(hvib)

    return Hvib



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
    dE = comn.energy_gaps(Hvib)

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


def run_LZ(params):
    """
    Main function to run the SH calculations based on the Landau-Zener hopping
    probabilities, all within the NBRA. The probabilities are implemented according to
    the Belyaev-Lebedev work.


    Required parameter keys:

    ===== Data description =====

    params["norbitals"]      [int] - how many lines/columns in the file
    params["active_space"]   [list of ints] - which orbitals we care about (indexing starts with 0)
    params["Hvib_re_prefix"] [string] - prefixes of the files with real part of the Hamiltonian
    params["Hvib_im_prefix"] [string] - prefixes of the files with imaginary part of the Hamiltonian
    params["Hvib_re_suffix"] [string] - suffixes of the files with real part of the Hamiltonian
    params["Hvib_im_suffix"] [string] - suffixes of the files with imaginary part of the Hamiltonian
    params["nsteps"]         [int] - how many files to read


    ===== Modeling params ===== 

    params["dt"]             [double, a.u.] - time distance between the adjacent data points
    params["ntraj"]          [int] - how many stochastic trajectories to use in the ensemble
    params["istate"]         [int] - index of the starting state (within those used in the active_space - see above)
    params["outfile"]        [string] - the name of the file, where all the results will be printed out
    params["T"]              [double, K] - temperature of the simulation

    """


    critical_params = [  ]
    default_params = { "T":300.0, "ntraj":1000, "nsteps":1,"istate":0, 
                       "sh_method":1, "decoherence_method":0, "dt":41.0,
                       "outfile":"_out.txt" }
    comn.check_input(params, default_params, critical_params)


    
    rnd = Random()

    #============ Read the data ===============
    Hvib = get_data(params)
    nstates= Hvib[0].num_of_cols

    #===== Precompute hopping probabilities ===
    P = Belyaev_Lebedev(Hvib, params["dt"])

    #======= Initialization ===================
    states = []
    for traj in xrange(params["ntraj"]):
        states.append(params["istate"])
    pops = tsh.compute_sh_statistics(nstates, states)


    #============== Dynamics ==================
    f = open(params["outfile"], "w");  f.close()

    for n in xrange(params["nsteps"]):

      comn.printout(n*params["dt"], pops, Hvib[n], params["outfile"])
      pops = tsh.compute_sh_statistics(nstates, states)
           

      for traj in xrange(params["ntraj"]):
          # Proposed hop:
          ksi = rnd.uniform(0.0, 1.0)
          st_new = tsh.hop_py(states[traj], P[n].T(), ksi)  

          # Accept the proposed hop with the Boltzmann probability
          ksi1 = rnd.uniform(0.0, 1.0)
          de = (Hvib[n].get(st_new,st_new) - Hvib[n].get(states[traj], states[traj])).real
          if de>0.0:
              if ksi1<prob.Boltz_quant_prob([0.0, de], params["T"])[1]:         
                  states[traj] = st_new
          else:
              states[traj] = st_new

