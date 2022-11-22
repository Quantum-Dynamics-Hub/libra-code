#*********************************************************************************                     
#* Copyright (C) 2022 Matthew Dutra and Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************

"""
..module:: initialize
  :platform: Unix, Windows
  :synopsis: This module contains functions for initial basis placement.

..moduleauthor :: Matthew Dutra, Alexey Akimov
"""

import sys
import os
import numpy as np

from liblibra_core import *
import util.libutil as comn

def initialize(_params):
    """Places the basis functions on all surfaces according to the `init_placement' parameter and
       initializes their qpas MATRIX objects.

    Args:
        _params (dict): Dictionary containing simulation parameters.

          * **_params[`init_placement`]** (float) : a parameter controlling how the initial placement
              of the basis functions is handled; 
          
            - 0: grid [ default ]
            - 1: gaussian distributed


    Also see:
        Parameters of the `grid` and `gaussian` functions


    Returns:
        (ntraj, Q, P, A, S, active_states): where:

        ntraj (int): The total number of trajectories across all surfaces. 
        Q (MATRIX(ndof, ntraj) ): coordinates of the GBFs
        P (MATRIX(ndof, ntraj) ): momenta of the GBFs
        A (MATRIX(ndof, ntraj) ): widths of the GBFs
        S (MATRIX(ndof, ntraj) ): total phases of the GBFs
        states (list of `ntraj` ints): quantum states for all trajectories

    """

    params = dict(_params)

    critical_params = [ ]
    default_params = { "init_placement":0 }
    comn.check_input(params, default_params, critical_params)

    init_placement = params['init_placement']

    if init_placement == 0:
        ntraj, Q, P, A, S, states = grid(params)

    elif init_placement == 1:
        ntraj, Q, P, A, S, states = gaussian(params)

    else:
        print("Unrecognized option in basis initialization! Available options: 0 = grid, 1 = gaussian.\n Exiting...\n")
        sys.exit(0)

    return ntraj, Q, P, A, S, states



def grid(_params):
    """Returns the initial basis parameters Q, P, A, S as ndof-by-ntraj matrices, 
       based on the input contained in the dict *dyn_params*. The placement is evenly spaced 
       across the domain where the initial wavefunction has a magnitude greater than *rho_cut*, and 
       each surface has the same number of trajectories. If *rho_cut* is negative, a user-defined
       set of boundary conditions specified by *grid_min* and *grid_max* is used instead.

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **_params[`states`]** (int) : states onto which to put trajectories [ default: [0, 1] ]

          * **_params[`grid_dims`]** (list of `ndof` floats) : the number of TBFs to be placed
              along each DOF for each surface. For grid, the list contains *ndof* elements, with each specifying the number
              of basis functions to be placed along each grid dimension. [ default: [5] ]

          * **_params[`rho_cut`]** (float) : cutoff parameter for basis placement, based on the wavefunction
              density [ default: 1e-12 ]

          * **_params[`alp_scl`]** (list of floats) : scaling parameter for the basis function width,
              based on the wavefunction width *wfc_a0* [ default: 8.0 ]

          * **_params[`wfc_q0`]** (list of floats) : list of coordinates (length *ndof*) for the initial 
              wavefunction position [ default: 0.0 ]

          * **_params[`wfc_p0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction momenta [ default: 0.0 ]

          * **_params[`wfc_a0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction widths [ default: 1.0 ]

          * **_params[`wfc_s0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction phase [ default: 0.0 ]
              
          * **_params[`grid_min`]** (list of floats) : list of minimum values for each domain when a user-defined cutoff is
              specified (supersedes the wavefunction-density based method)
              
          * **_params[`grid_max`]** (list of floats) : list of maximum values for each domain when a user-defined cutoff is
              specified (supersedes the wavefunction-density based method)


    Returns:
        (ntraj, Q, P, A, S, active_states): where:

        ntraj (int): The total number of trajectories across all surfaces. 
        Q (MATRIX(ndof, ntraj) ): coordinates of the GBFs
        P (MATRIX(ndof, ntraj) ): momenta of the GBFs
        A (MATRIX(ndof, ntraj) ): widths of the GBFs
        S (MATRIX(ndof, ntraj) ): total phases of the GBFs
        states (list of `ntraj` ints): quantum states for all trajectories
    """

    params = dict(_params)

    critical_params = [ ]
    default_params = { "grid_dims":[5], "rho_cut":1e-12, "alp_scl":[8.0], 
                       "wfc_q0":[0.0], "wfc_p0":[0.0], "wfc_a0":[1.0], "wfc_s0":[0.0],
                       "grid_min":[-1.0], "grid_max":[1.0]
                     }
    comn.check_input(params, default_params, critical_params)

    states = params['states']
    grid_dims = params['grid_dims']
    rho_cut = params['rho_cut']
    q0 = params['wfc_q0']
    p0 = params['wfc_p0']
    a0 = params['wfc_a0']
    alp_scl = params['alp_scl']
    s0 = params['wfc_s0']
    grid_min = params["grid_min"]
    grid_max = params["grid_max"]

    nstates = len(states)
    ndof = len(q0)

    ntraj_on_state = 1
    for i in range(len(grid_dims)):
        ntraj_on_state *= grid_dims[i]

    ntraj = ntraj_on_state*nstates


    qvals = MATRIX(ndof,ntraj)
    pvals = MATRIX(ndof,ntraj)
    avals = MATRIX(ndof,ntraj)
    svals = MATRIX(ndof,ntraj)
    surf_ids = []


    qlo,qhi = [], []

    if rho_cut > 0.0:
        for dof in range(ndof):
            dq = np.sqrt(-0.5/a0[dof]*np.log(rho_cut))
            xlow = q0[dof] - dq
            xhi  = q0[dof] + dq
            qlo.append(xlow)
            qhi.append(xhi)
    else:
        qlo = list(grid_min)
        qhi = list(grid_max)
    

    grid_bounds = np.mgrid[tuple(slice(qlo[dof],qhi[dof],complex(0,grid_dims[dof])) for dof in range(ndof))]

    elems = 1
    for i in grid_dims:
        elems *= i

    b = grid_bounds.flatten()
    qs = []
    for i in range(elems):
        index = i
        coords = []

        while index < len(b):
            coords.append(b[index])
            index += elems
        qs.append(coords)

    for dof in range(ndof):
        for j in range(ntraj_on_state):
            for n in range(nstates):
                qvals.set(dof,j+n*ntraj_on_state, qs[j][dof])
                pvals.set(dof,j+n*ntraj_on_state, p0[dof])
                avals.set(dof,j+n*ntraj_on_state, a0[dof]*alp_scl[dof])
                svals.set(dof,j+n*ntraj_on_state,0.0)

    for n in states:
        for j in range(ntraj_on_state):
            surf_ids.append(n)


    return ntraj, qvals, pvals, avals, svals, surf_ids


def gaussian(_params):
    """Returns the initial basis parameters Q, P, A, S as ndof-by-ntraj matrices,
       based on the input contained in the dict *dyn_params*. The placement is randomly chosen from a 
       Gaussian distribution centered about the wavepacket maximum with a standard deviation *rho_cut*, 
       and each surface has the same number of trajectories. The corresponding Gaussian-distributed 
       momenta are ordered so that the wavepacket spreads as x increases.

    Args:
        _params (dict): Dictionary containing simulation parameters.

          * **_params[`states`]** (list of int) : states onto which put trajectories [ default: [0, 1] ]

          * **_params[`grid_dims`]** (list of floats) : the total number of basis functions to be 
              placed on each surface. For Gaussian, the list has only one element, specifying the 
              number of basis functions per surface. Note that the total number of basis functions 
              will then be *nstates*-by-*prod(grid_dims)*

          * **_params[`rho_cut`]** (float) : cutoff parameter for basis placement, based on the 
              wavefunction density (strictly speaking, the standard deviation of the Gaussian sampling)
              [ default: 1e-12 ]

          * **_params[`alp_scl`]** (list of floats) : scaling parameter for the basis function width,
              based on the wavefunction width *wfc_a0* [ default: 8.0 ]

          * **_params[`wfc_q0`]** (list of floats) : list of coordinates (length *ndof*) for the initial
              wavefunction position [ default: 0.0 ]

          * **_params[`wfc_p0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction momenta [ default: 0.0 ]

          * **_params[`wfc_a0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction widths [ default: 1.0 ]

          * **_params[`wfc_s0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction phase [ default: 0.0 ]

    Returns:
        (ntraj, Q, P, A, S, active_states): where:

        ntraj (int): The total number of trajectories across all surfaces. 
        Q (MATRIX(ndof, ntraj) ): coordinates of the GBFs
        P (MATRIX(ndof, ntraj) ): momenta of the GBFs
        A (MATRIX(ndof, ntraj) ): widths of the GBFs
        S (MATRIX(ndof, ntraj) ): total phases of the GBFs
        states (list of `ntraj` ints): quantum states for all trajectories

    """

    params = dict(_params)

    critical_params = []
    default_params = { "grid_dims":[5], "rho_cut":1e-12, 
                       "wfc_q0":[0.0], "wfc_p0":[0.0], "wfc_a0":[1.0], "alp_scl":[1.0], "wfc_s0":[0.0]
                     }
    comn.check_input(params, default_params, critical_params)



    states = params["states"]
    nstates = len(states)
    grid_dims = params['grid_dims']
    rho_cut = params['rho_cut']
    q0 = params['wfc_q0']
    p0 = params['wfc_p0']
    a0 = params['wfc_a0']
    alp_scl = params['alp_scl']
    s0 = params['wfc_s0']

    ndof = len(q0)

    ntraj_on_state = 1
    for i in range(len(grid_dims)):
        ntraj_on_state *= grid_dims[i]

    ntraj = ntraj_on_state*nstates
    qvals=MATRIX(ndof,ntraj)
    pvals=MATRIX(ndof,ntraj)
    avals=MATRIX(ndof,ntraj)
    svals=MATRIX(ndof,ntraj)

    surf_ids = []
    ntraj = grid_dims[0]*nstates

    q_gaus, p_gaus = [],[]
    for dof in range(ndof):
        q_gaus.append(np.sort(np.random.normal(q0[dof], rho_cut, ntraj)))
        p_gaus.append(np.sort(np.random.normal(p0[dof], rho_cut, ntraj)))

    for dof in range(ndof):
        for traj in range(ntraj):
            qvals.set(dof,traj, q_gaus[dof][traj])
            pvals.set(dof,traj, p_gaus[dof][traj])
            avals.set(dof,traj, a0[dof]*alp_scl[dof])
            svals.set(dof,traj, 0.0)

    for n in states:
        for j in range(ntraj_on_state):
            surf_ids.append(n)

    return ntraj, qvals, pvals, avals, svals, surf_ids


def coeffs(aQ, aP, aA, aS, astate, bQ, bP, bA, bS, bstate):
    """Returns the projection vector C of the initial wavefunction |psi> onto the basis of GBFs |G>={|G_i>}: 

       |psi> = sum_i {  |G_i> * C_i } = |G> * C, so:

       C = S^{-1} * <G|psi>

       The initial (active, so "a") wavefunction is given by a multidimensional Gaussian with the parameters aQ, aP, aA, aS, astate
       The GBFs (basis, so "b") are parameterized by bQ, bP, bA, bS, bstate

    Args:

        aQ (MATRIX(ndof, 1) ): coordinate of the initial wavefunction, |psi>
        aP (MATRIX(ndof, 1) ): momenta of the initial wavefunction, |psi>
        aA (MATRIX(ndof, 1) ): widths of the initial wavefunction, |psi>
        aS (MATRIX(ndof, 1) ): phases of the initial wavefunction, |psi>
        astate (list of `ntraj` ints): quantum states the initial wavefunction, |psi>
        bQ (MATRIX(ndof, ntraj) ): coordinates of the GBFs
        bP (MATRIX(ndof, ntraj) ): momenta of the GBFs
        bA (MATRIX(ndof, ntraj) ): widths of the GBFs
        bS (MATRIX(ndof, ntraj) ): phases of the GBFs
        bstate (list of `ntraj` ints): quantum states for all trajectories

    Returns:

        (CMATRIX): Projection vector C for the initial wavefunction onto the initial Gaussian basis.

    """

    ndof = bQ.num_of_rows
    ntraj = bQ.num_of_cols
    nstates = len(set(bstate))

    bA_half = 0.5*bA
    
    # Not sure why we were considering full A for the initial wavefunction
    G_psi = gwp_overlap_matrix(bQ, bP, bS, bA_half, Py2Cpp_int(bstate), aQ, aP, aS, aA, Py2Cpp_int(astate))

    GG = gwp_overlap_matrix(bQ, bP, bS, bA_half, Py2Cpp_int(bstate), bQ, bP, bS, bA_half, Py2Cpp_int(bstate))

    tmp = FullPivLU_rank_invertible(GG)
    if tmp[1]==0:
        print("GBF overlap matrix is not invertible.\n Exiting...\n")
        sys.exit(0)

    invGG = CMATRIX(ntraj, ntraj)
    FullPivLU_inverse(GG, invGG)

    coeff = invGG * G_psi

    return coeff      


