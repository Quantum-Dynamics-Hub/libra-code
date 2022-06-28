#*********************************************************************************
#* Copyright (C) 2021-2022 Matthew Dutra, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#***********************************************************************************
"""
..module:: qtag_prop
  :platform: Unix, Windows
  :synopsis: This module contains functions for basis propagation in different multi-surface schemes (mss).

..moduleauthors :: Matthew Dutra, Alexey Akimov
"""

import sys
import os
from liblibra_core import *
from libra_py import data_outs

import numpy as np
from . import qtag_calc
from . import qtag_mom

import util.libutil as comn

def propagate(dyn_params, qpas, coeff, surf_pops):
    """Makes the trajectories on surfaces with low populations (where quantum momentum
       is ill-defined) move as the trajectories on the surface with the highest population

    Args:

        dyn_params (dict): Dictionary containing simulation parameters.
 
          * **dyn_params[`decpl_den`]** (float) : a parameter controlling independence of trajectories. If the population
              on a given state is larger that this parameter, the trajectories evolve according to their own
              quantum momentum. If the population is less than this number, the trajectories evolve according
              to the quantum momentum for the most populated state.  [ default: 0.25 ]

        qpas (list): List of {q,p,a,s} MATRIX objects.

        coeff (CMATRIX(ntraj x 1)): The complex coefficient matrix for the TBF on all surfaces.

        surf_pops (list): List of surface populations.

    Returns:
        qpas_new (list): List of updated {q,p,a,s} MATRIX objects for active surface

        btot (CMATRIX(ntraj x 1)): The complex projection vector for the TBF on all surfaces.
    """

    params = dict(dyn_params)

    critical_params = [  ]
    default_params = {"decpl_den":0.25}
    comn.check_input(params, default_params, critical_params)

    dt = params["dt"]
    decpl = params["decpl_den"] 
    iM = params["iM"]  # MATRIX(ndof, 1)
    states = params["states"]

    q_update_method = params["q_update_method"]  # 0 - frozen, 1 - move
    p_update_method = params["p_update_method"]  # 0 - frozen, 1 - move
    a_update_method = params["a_update_method"]  # 0 - frozen, 1 - move
    s_update_method = params["s_update_method"]  # 0 - frozen, 1 - move

    q_sync_method = params["q_sync_method"] # 0 - use the current value, 1 - use the value from the most populated surface
    p_sync_method = params["p_sync_method"] # 0 - use the current value, 1 - use the value from the most populated surface
    a_sync_method = params["a_sync_method"] # 0 - use the current value, 1 - use the value from the most populated surface
    s_sync_method = params["s_sync_method"] # 0 - use the current value, 1 - use the value from the most populated surface

    # The original set
    q_old = MATRIX(qpas[0])
    p_old = MATRIX(qpas[1])
    a_old = MATRIX(qpas[2])
    s_old = MATRIX(qpas[3])    
    surf_ids = qpas[4]

    ndof = q_old.num_of_rows
    ntraj = q_old.num_of_cols
    nstates = len(states)

    x_dofs = list(range(ndof))
    invM = MATRIX(ndof, ndof)
    for i in x_dofs:
        invM.set(i, i,  iM.get(i, 0))

    # The new set (propagated)
    q_new = MATRIX(q_old)
    p_new = MATRIX(p_old)
    a_new = MATRIX(a_old)
    s_new = MATRIX(s_old)

#    int ii = 0
#    unsorted_pairs = []
#    for i in states:
#        unsorted_pairs.append([n, surf_pops[i]])
#        ii += 1

#    sorted_pairs = merge_sort(unsorted_pairs) 

#    sorted_states = [0]*nstates
#    for i in range(nstates):
#        sorted_states[nstates-1-i] = sorted_pairs[i][0]

    sorted_pops = sorted(surf_pops, reverse = True)
    sorted_states = []
    for n in range(nstates):
        indx = surf_pops.index(sorted_pops[n])
        sorted_states.append(states[indx])

   
    # The properties for the trajectories on the most populated surface - reference for synchronizing
    q_new_on_surf_ref = None
    p_new_on_surf_ref = None
    a_new_on_surf_ref = None
    s_new_on_surf_ref = None

#    iref_pop_state = True
    for nindex, n in enumerate(sorted_states): # over all states, but starting with the most populated one

        traj_on_surf = [index for index, traj_id in enumerate(surf_ids) if traj_id == n]
        ntraj_on_surf = len(traj_on_surf)  
        surf_pop_ref = sorted_pops[nindex]

        q_on_surf = MATRIX(ndof, ntraj_on_surf)
        p_on_surf = MATRIX(ndof, ntraj_on_surf)
        a_on_surf = MATRIX(ndof, ntraj_on_surf)
        s_on_surf = MATRIX(ndof, ntraj_on_surf)

        pop_submatrix(q_old, q_on_surf, x_dofs, traj_on_surf)
        pop_submatrix(p_old, p_on_surf, x_dofs, traj_on_surf)
        pop_submatrix(a_old, a_on_surf, x_dofs, traj_on_surf)
        pop_submatrix(s_old, s_on_surf, x_dofs, traj_on_surf)

#Forcing sync'ed trajectories for q_sync_method = 2, regardless of pop
#        if q_sync_method == 2:
#            if n>0:
#                surf_pop_ref = 0.0

        # "Independent" evolution of trajectories
        if surf_pop_ref > decpl:
            coeff_on_surf = CMATRIX(ntraj_on_surf,1)  # coefficients for the TBFs on the surface n
            pop_submatrix(coeff, coeff_on_surf, traj_on_surf, [0])

            qpas_on_surf = [q_on_surf, p_on_surf, a_on_surf, s_on_surf] # variables for the TBFs on the surface n
            mom, r, gmom, gr = qtag_mom.momentum(dyn_params, qpas_on_surf, coeff_on_surf)

            # mom - MATRIX(ndof, ntraj_on_surf) - quantum momentum for q -  Im( nabla_{\alp} \psi / \psi)
            # r - MATRIX(ndof, ntraj_on_surf) - quantum momentum for s   -  Re( nabla_{\alp} \psi / \psi)
            # gmom - MATRIX(ndof, ntraj_on_surf) - gradient of the fit of mom ~ nabla_{\alp} Im( nubla_{\alp} \psi / \psi) - to update alphas
            # gr - MATRIX(ndof, ntraj_on_surf) - gradient of the fit of r ~ nubla_{\alp} Re( nubla_{\alp} \psi / \psi)  - to update alphas

            # mom_tmp = qtag_momentum(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff);
            # mom = mom_tmp.imag() 
            # r = mom_tmp.real()

            if q_update_method==0:
                q_new_on_surf = MATRIX(q_on_surf)
            elif q_update_method==1:
                q_new_on_surf = q_on_surf + invM * mom * dt 
            else:
                print(F"q_update_method = {q_update_method} is not implemented. Exiting...\n")
                sys.exit(0)

            if p_update_method==0:
                p_new_on_surf = MATRIX(p_on_surf)
            elif p_update_method==1:
                p_new_on_surf = MATRIX(mom)
            else:
                print(F"p_update_method = {p_update_method} is not implemented. Exiting...\n")
                sys.exit(0)


            if a_update_method==0:
                a_new_on_surf = MATRIX(a_on_surf)
            elif a_update_method==1:
                a_tmp = MATRIX(ndof, ntraj_on_surf)
                a_tmp.dot_product(a_on_surf, gmom)
                a_new_on_surf = a_on_surf - 2.0*dt * invM * a_tmp
            else:
                print(F"a_update_method = {a_update_method} is not implemented. Exiting...\n")
                sys.exit(0)

            s_new_on_surf = MATRIX(s_on_surf)

            # Put the surface-specific variables in the global variables
            push_submatrix(q_new, q_new_on_surf, x_dofs, traj_on_surf)
            push_submatrix(p_new, p_new_on_surf, x_dofs, traj_on_surf)
            push_submatrix(a_new, a_new_on_surf, x_dofs, traj_on_surf)
            push_submatrix(s_new, s_new_on_surf, x_dofs, traj_on_surf)

             # Save it for later - for the "dependent" trajectories
            if nindex == 0:
                q_new_on_surf_ref = MATRIX(q_new_on_surf)
                p_new_on_surf_ref = MATRIX(p_new_on_surf)
                a_new_on_surf_ref = MATRIX(a_new_on_surf)
                s_new_on_surf_ref = MATRIX(s_new_on_surf)

        else:
            # For this to work, we need that all surfaces have equal number of trajectories
            if q_sync_method >= 1:
                push_submatrix(q_new, q_new_on_surf_ref, x_dofs, traj_on_surf)
            if p_sync_method >= 1:
                push_submatrix(p_new, p_new_on_surf_ref, x_dofs, traj_on_surf)
            if a_sync_method >= 1:
                push_submatrix(a_new, a_new_on_surf_ref, x_dofs, traj_on_surf)
            if s_sync_method >= 1:
                push_submatrix(s_new, s_new_on_surf_ref, x_dofs, traj_on_surf)


    qpas_new = [q_new, p_new, a_new, s_new, surf_ids]

#    ov_no = qtag_calc.new_old_overlap(ndof, ntraj, states, qpas, qpas_new)
    st = qtag_calc.time_overlap(ndof, ntraj, states, qpas_new, qpas)
    btot = st*coeff

    return qpas_new, btot


