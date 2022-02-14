"""
..module:: qtag_init
  :platform: Unix, Windows
  :synopsis: This module contains functions for initial basis placement.

..moduleauthors :: Matthew Dutra
"""

import sys
import os
from liblibra_core import *

import numpy as np

def grid(ndof,nstates,params):
    """Returns the initial basis parameters {q,p,a,s} as a list of  ndof-by-ntraj matrices *qpas*, 
       based on the input contained in the dicts *traj0* and *wf0*. The placement is evenly spaced 
       across the domain where the initial wavefunction has a magnitude greater than *rcut*, and 
       each surface has the same number of trajectories.

    Args:
        ndof (integer): The number of degrees of freedom.

        nstates (integer): The total number of states.

        traj0 (dictionary): Dictionary containing initialization parameters and keywords.

        wf0 (dictionary): Dictionary containing initial wavepacket conditions.

    Returns:
        ntraj (integer): The number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects storing trajectory information.
    """

    grid_dims=params['grid_dims']

    ntraj_on_state = 1
    for i in range(len(grid_dims)):
        ntraj_on_state *= grid_dims[i]

    ntraj = ntraj_on_state*nstates

    qvals=MATRIX(ndof,ntraj)
    pvals=MATRIX(ndof,ntraj)
    avals=MATRIX(ndof,ntraj)
    svals=MATRIX(ndof,ntraj)

    rcut=params['rho_cut']
    a0=params['a0']

    surf_ids = []
    qlo,qhi = [], []

    for dof in range(ndof):
        xlow = wf0['q'][dof]-np.sqrt(-0.5/wf0['a'][dof]*np.log(rcut))
        xhi = params['basis_qtype'][dof]+np.sqrt(-0.5/params['wfa0'][dof]*np.log(rcut))
        qlo.append(xlow)
        qhi.append(xhi)

    grid = np.mgrid[tuple(slice(qlo[dof],qhi[dof],complex(0,grid_dims[dof])) for dof in range(ndof))]

    elems = 1
    for i in grid_dims:
        elems *= i

    b = grid.flatten()
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
                qvals.set(dof,j+n*ntraj_on_state,qs[j][dof])
                pvals.set(dof,j+n*ntraj_on_state,params['wfp0'][dof])
                avals.set(dof,j+n*ntraj_on_state,params['wfa0'][dof]*a0[dof])
                svals.set(dof,j+n*ntraj_on_state,0.0)

    for n in range(nstates):
        for j in range(ntraj_on_state):
            surf_ids.append(n)

    qpas=[qvals,pvals,avals,svals,surf_ids]
    return(ntraj,qpas)

def gaussian(ndof,ntraj,traj0,wf0):
    """Returns the initial basis parameters {q,p,a,s} as an *ntraj*-by-4 matrix *qpas*, based on the input contained in the dicts *traj0* and *wf0*. The placement is randomly chosen from a Gaussian distribution centered about the wavepacket maximum with a standard deviation *rho*. The corresponding Gaussian-distributed momenta are ordered so that the wavepacket spreads as x increases.

    Args:
        ndof (integer): The number of degrees of freedom.

        ntraj (integer): The number of trajectories per surface.

        traj0 (dictionary): Dictionary containing initialization parameters and keywords.

        wf0 (dictionary): Dictionary containing initial wavepacket conditions.

    Returns:
        qpas (list): List of {q,p,a,s} MATRIX objects.
    """

    qvals=MATRIX(ndof,ntraj)
    pvals=MATRIX(ndof,ntraj)
    avals=MATRIX(ndof,ntraj)
    svals=MATRIX(ndof,ntraj)

    rho=traj0['rho']
    a0=traj0['a0']

    q_gaus=np.sort(np.random.normal(wf0['q'],rho,ntraj))
    p_gaus=np.sort(np.random.normal(wf0['p'],rho,ntraj))

    for i in range(ndof):
        for j in range(ntraj):
            qvals.set(i,j,q_gaus[i])
            pvals.set(i,j,p_gaus[i])
            avals.set(i,j,wf0['a']*a0)
            svals.set(i,j,0.0)

    qpas=[qvals,pvals,avals,svals]
    return(qpas)

def coeffs(params,qpas,active_state):
    """Returns the projection vector *b* of the initial wavefunction with parameters stored in the dict *wf0* onto the basis defined by *qpas*. This function assumes the wavefunction is located entirely on *nsurf*=1 initially.

    Args:
        wf0 (dictionary): Dictionary containing initial wavepacket conditions.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        active_state (integer): Number specifying the active state (i.e. the state w/ the initial wavepacket). 0 = ground, 1 = first excited, ...

    Returns:

        b (CMATRIX): Projection vector for the initial wavefunction onto the initial Gaussian basis.
    """

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
        q2.set(dof,0,params['wfq0'][dof])
        p2.set(dof,0,params['wfp0'][dof])
        a2.set(dof,0,params['wfa0'][dof])
        s2.set(dof,0,params['wfs0'][dof])

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

