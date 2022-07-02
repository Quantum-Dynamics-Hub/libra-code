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
..module:: qtag_init
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
    """Returns the initial basis parameters {q,p,a,s} as a list of  ndof-by-ntraj matrices *qpas*, 
       based on the input contained in the dict *dyn_params*. The placement is evenly spaced 
       across the domain where the initial wavefunction has a magnitude greater than *rho_cut*, and 
       each surface has the same number of trajectories.

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
                       "wfc_q0":[0.0], "wfc_p0":[0.0], "wfc_a0":[1.0], "wfc_s0":[0.0]
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

    for n in states:
        for j in range(ntraj_on_state):
            surf_ids.append(n)


    return ntraj, qvals, pvals, avals, svals, surf_ids


def gaussian(_params):
    """Returns the initial basis parameters {q,p,a,s} as a list of  ndof-by-ntraj matrices *qpas*,
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


#def restart(dyn_params):

#    prefix = dyn_params['prefix']
#    states = dyn_params['states']
#    grid_dims = dyn_params['grid_dims']
    
#    nstates = len(states)

#    ntraj_on_state = 1
#    for i in range(len(grid_dims)):
#        ntraj_on_state *= grid_dims[i]

#    ntraj = ntraj_on_state*nstates

#    qfile = f"{prefix}/q.txt"
#    try:
#        os.path.exists(qfile)
#    except:
#        sys.exit("ERROR in qtag_init.restart() -- necessary file 'q.txt' is missing from "+prefix+" directory!"

#    pfile = f"{prefix}/p.txt"
#    try: 
#        os.path.exists(pfile)
#    except:
#        sys.exit("ERROR in qtag_init.restart() -- necessary file 'p.txt' is missing from "+prefix+" directory!"

#    afile = f"{prefix}/a.txt"
#    try: 
#        os.path.exists(afile)
#    except:
#        sys.exit("ERROR in qtag_init.restart() -- necessary file 'a.txt' is missing from "+prefix+" directory!"

#    sfile = f"{prefix}/s.txt"
#    try: 
#        os.path.exists(sfile)
#    except:
#        sys.exit("ERROR in qtag_init.restart() -- necessary file 's.txt' is missing from "+prefix+" directory!"

#    cfile = f"{prefix}/coeffs.txt"
#    try: 
#        os.path.exists(cfile)
#    except:
#        sys.exit("ERROR in qtag_init.restart() -- necessary file 'coeffs.txt' is missing from "+prefix+" directory!"

#    qdata = data_read.get_data_from_file2(qfile, [traj for traj in range(ntraj)])
#    pdata = data_read.get_data_from_file2(pfile, [traj for traj in range(ntraj)])
#    adata = data_read.get_data_from_file2(afile, [traj for traj in range(ntraj)])
#    sdata = data_read.get_data_from_file2(sfile, [traj for traj in range(ntraj)])
#    cdata = data_read.get_data_from_file2(cfile, [traj for traj in range(ntraj)])

#    return ntraj,qpas


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

    # qvals - basis, q2 - initial wfc
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
    c0 = CMATRIX(ntraj,1)

    for dof in range(ndof):
        q2.set(dof,0, q0[dof])
        p2.set(dof,0, p0[dof])
        a2.set(dof,0, a0[dof])
        s2.set(dof,0, s0[dof])

    for i in range(ntraj):
        b.set(i,complex(0.0,0.0))
        c0.set(i,complex(0.0,0.0))

    traj_active = [index for index, traj_id in enumerate(surf_ids) if traj_id == active_state]
    ntraj_active = len(traj_active)

    qvals_active = MATRIX(ndof,ntraj_active)
    pvals_active = MATRIX(ndof,ntraj_active)
    avals_active = MATRIX(ndof,ntraj_active)
    svals_active = MATRIX(ndof,ntraj_active)
    bvals_active = CMATRIX(ntraj_active,1)
    ctemp = CMATRIX(ntraj_active,1)

    basis_ovlp = CMATRIX(ntraj_active,ntraj_active)

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

    for i in range(ntraj_active):
        qi = qvals_active.col(i)
        pi = pvals_active.col(i)
        ai = avals_active.col(i)
        si = svals_active.col(i)

        for j in range(ntraj_active):
            qj = qvals_active.col(j)
            pj = pvals_active.col(j)
            aj = avals_active.col(j)
            sj = svals_active.col(j)
            basis_ovlp.set(i,j,gwp_overlap(qi,pi,si,ai/2,qj,pj,sj,aj/2))

#    basis_ovlp = gwp_overlap_matrix(qvals_active, pvals_active, svals_active, 0.5*avals_active, 
#                                    qvals_active, pvals_active, svals_active, 0.5*avals_active)

    pop_submatrix(b,bvals_active,traj_active,[0])
    ctemp = basis_ovlp * bvals_active

    ii = 0
    for i in traj_active:
        c0.set(i,0,ctemp.get(ii))
        ii += 1

    return (b, c0)
    """


