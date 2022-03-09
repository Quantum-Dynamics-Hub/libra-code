"""
..module:: qtag_init
  :platform: Unix, Windows
  :synopsis: This module contains functions for initial basis placement.

..moduleauthors :: Matthew Dutra
"""

import sys
import os
import numpy as np

from liblibra_core import *
import util.libutil as comn


def initialize(dyn_params):
    """Places the basis functions on all surfaces according to the `init_placement' parameter and
       initializes their qpas MATRIX objects.

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`init_placement`]** (float) : a parameter controlling how the initial placement
              of the basis functions is handled; 0 = grid, 1 = gaussian distributed [ default: 0 ]

    Returns:
        ntraj (int): The total number of trajectories across all surfaces.

        qpas (list): List of {q,p,a,s} MATRIX objects.
    """

    params = dict(dyn_params)

    critical_params = [ ]
    default_params = { "init_placement":0 }
    comn.check_input(params, default_params, critical_params)

    init_placement = params['init_placement']

    if init_placement == 0:
        ntraj, qpas = grid(dyn_params)

    elif init_placement == 1:
        ntraj, qpas = gaussian(dyn_params)

    else:
        print("Unrecognized option in basis initialization! 0 = grid, 1 = gaussian.\n Exiting...\n")
        sys.exit(0)

    return ntraj, qpas


def grid(dyn_params):
    """Returns the initial basis parameters {q,p,a,s} as a list of  ndof-by-ntraj matrices *qpas*, 
       based on the input contained in the dict *dyn_params*. The placement is evenly spaced 
       across the domain where the initial wavefunction has a magnitude greater than *rho_cut*, and 
       each surface has the same number of trajectories.

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`nstates`]** (int) : the number of states [ default: 2 ]

          * **dyn_params[`grid_dims`]** (list of floats) : the total number of basis functions to be placed
              on each surface. For grid, the list contains *ndof* elements, with each specifying the number
              of basis functions to be placed along each grid dimension. [ default: 5 ]

          * **dyn_params[`rho_cut`]** (float) : cutoff parameter for basis placement, based on the wavefunction
              density [ default: 1e-12 ]

          * **dyn_params[`alp_scl`]** (list of floats) : scaling parameter for the basis function width,
              based on the wavefunction width *wfc_a0* [ default: 8.0 ]

          * **dyn_params[`wfc_q0`]** (list of floats) : list of coordinates (length *ndof*) for the initial 
              wavefunction position [ default: 0.0 ]

          * **dyn_params[`wfc_p0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction momenta [ default: 0.0 ]

          * **dyn_params[`wfc_a0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction widths [ default: 1.0 ]

          * **dyn_params[`wfc_s0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction phase [ default: 0.0 ]

    Returns:
        ntraj (integer): The number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects storing trajectory information.
    """

    params = dict(dyn_params)

    critical_params = [ "nstates" ]
    default_params = { "grid_dims":[5], "rho_cut":1e-12, "alp_scl":[8.0], 
                       "wfc_q0":[0.0], "wfc_p0":[0.0], "wfc_a0":[1.0], "wfc_s0":[0.0]
                     }
    comn.check_input(params, default_params, critical_params)

    nstates = params["nstates"]
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
    qlo,qhi = [], []

    for dof in range(ndof):
        dq = np.sqrt(-0.5/a0[dof]*np.log(rho_cut))
        xlow = q0[dof] - dq
        xhi  = q0[dof] + dq
        qlo.append(xlow)
        qhi.append(xhi)

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

    for n in range(nstates):
        for j in range(ntraj_on_state):
            surf_ids.append(n)

    qpas=[qvals,pvals,avals,svals,surf_ids]

    return ntraj, qpas


def gaussian(dyn_params):
    """Returns the initial basis parameters {q,p,a,s} as a list of  ndof-by-ntraj matrices *qpas*,
       based on the input contained in the dict *dyn_params*. The placement is randomly chosen from a 
       Gaussian distribution centered about the wavepacket maximum with a standard deviation *rho_cut*, 
       and each surface has the same number of trajectories. The corresponding Gaussian-distributed 
       momenta are ordered so that the wavepacket spreads as x increases.

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`nstates`]** (int) : the number of states [ default: 2 ]

          * **dyn_params[`grid_dims`]** (list of floats) : the total number of basis functions to be 
              placed on each surface. For Gaussian, the list has only one element, specifying the 
              number of basis functions per surface. Note that the total number of basis functions 
              will then be *nstates*-by-*prod(grid_dims)*

          * **dyn_params[`rho_cut`]** (float) : cutoff parameter for basis placement, based on the 
              wavefunction density (strictly speaking, the standard deviation of the Gaussian sampling)
              [ default: 1e-12 ]

          * **dyn_params[`alp_scl`]** (list of floats) : scaling parameter for the basis function width,
              based on the wavefunction width *wfc_a0* [ default: 8.0 ]

          * **dyn_params[`wfc_q0`]** (list of floats) : list of coordinates (length *ndof*) for the initial
              wavefunction position [ default: 0.0 ]

          * **dyn_params[`wfc_p0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction momenta [ default: 0.0 ]

          * **dyn_params[`wfc_a0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction widths [ default: 1.0 ]

          * **dyn_params[`wfc_s0`]** (list of floats) : list of values (length *ndof*) for the initial
              wavefunction phase [ default: 0.0 ]

    Returns:
        ntraj (integer): The number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects storing trajectory information.
    """

    params = dict(dyn_params)

    critical_params = [ "nstates" ]
    default_params = { "grid_dims":[5], "rho_cut":1e-12, 
                       "wfc_q0":[0.0], "wfc_p0":[0.0], "wfc_a0":[1.0], "alp_scl":[1.0], "wfc_s0":[0.0]
                     }
    comn.check_input(params, default_params, critical_params)



    nstates = params["nstates"]
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

    for n in range(nstates):
        for j in range(ntraj_on_state):
            surf_ids.append(n)

    qpas=[qvals,pvals,avals,svals,surf_ids]
    return ntraj,qpas



def coeffs(dyn_params, qpas, active_state):
    """Returns the projection vector *b* of the initial wavefunction with parameters stored in the dict 
       *dyn_params* onto the basis defined by *qpas*. This function assumes the wavefunction is initialized 
       on the state specified by *active_state*.

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        active_state (integer): Number specifying the active state (i.e. the state w/ the initial wavepacket).
        0 = ground, 1 = first excited, ...

    Returns:

        b (CMATRIX): Projection vector for the initial wavefunction onto the initial Gaussian basis.
    """


    params = dict(dyn_params)

    critical_params = [  ]
    default_params = { "wfc_q0":[0.0], "wfc_p0":[0.0], "wfc_a0":[1.0], "alp_scl":[1.0], "wfc_s0":[0.0]
                     }
    comn.check_input(params, default_params, critical_params)


#    nstates = params["nstates"]
    q0 = params['wfc_q0']
    p0 = params['wfc_p0']
    a0 = params['wfc_a0']
    s0 = params['wfc_s0']

    ndof = len(q0)



    qvals = qpas[0]
    pvals = qpas[1]
    avals = qpas[2]
    svals = qpas[3]
    surf_ids = qpas[4]

    ndof = qvals.num_of_rows
    ntraj = qvals.num_of_cols
    nstates = len(set(surf_ids))

    q2 = MATRIX(ndof,1)
    p2 = MATRIX(ndof,1)
    a2 = MATRIX(ndof,1)
    s2 = MATRIX(ndof,1)
    b = CMATRIX(ntraj,1)

    for dof in range(ndof):
        q2.set(dof,0, q0[dof])
        p2.set(dof,0, p0[dof])
        a2.set(dof,0, a0[dof])
        s2.set(dof,0, s0[dof])

    for i in range(ntraj):
        b.set(i,complex(0.0,0.0))

    traj_active = [index for index, traj_id in enumerate(surf_ids) if traj_id == active_state]
    ntraj_active = len(traj_active)

    qvals_active = MATRIX(ndof,ntraj_active)
    pvals_active = MATRIX(ndof,ntraj_active)
    avals_active = MATRIX(ndof,ntraj_active)
    svals_active = MATRIX(ndof,ntraj_active)

    pop_submatrix(qvals,qvals_active,[dof for dof in range(ndof)],traj_active)
    pop_submatrix(pvals,pvals_active,[dof for dof in range(ndof)],traj_active)
    pop_submatrix(avals,avals_active,[dof for dof in range(ndof)],traj_active)
    pop_submatrix(svals,svals_active,[dof for dof in range(ndof)],traj_active)

    ii = 0
    for i in traj_active:

        q1 = qvals_active.col(ii)
        p1 = pvals_active.col(ii)
        a1 = avals_active.col(ii)
        s1 = svals_active.col(ii)

        b.set(i,0,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2))
        ii += 1

    return(b)

