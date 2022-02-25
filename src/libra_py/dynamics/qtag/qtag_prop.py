"""
..module:: qtag_prop
  :platform: Unix, Windows
  :synopsis: This module contains functions for basis propagation in different multi-surface schemes (mss).

..moduleauthors :: Matthew Dutra
"""

import sys
import os
from liblibra_core import *
from libra_py import data_outs

import numpy as np
from . import qtag_calc
from . import qtag_mom


def propagate(params, qpas, coeff, surf_pops):
    """Makes the trajectories on surfaces with low populations (where quantum momentum
       is ill-defined) to move as the trajectories on the surface with the highest population

       Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), 
       as well as their corresponding projection vectors *b1* and *b2*, where the motion of both sets of 
       functions are synced to the lower energetic surface while the density on the upper surface is less than 
       a threshold value specified by the *decpl* parameter. Also necessary are the functions for calculating 
       momentum (*mom_calc*) and basis updates (*props*).

    Args:

        params (dict): Dictionary containing simulation parameters.
 
          * **params[`decpl_den`]** (float) : a parameter controlling independence of trajectories. If the population
              on a given state is larger that this parameter, the trajectories evolve according to their own
              quantum momentum. If the population is less than this number, the trajectories evolve according
              to the quantum momentum for the most populated state.  [ default: 0.3 ]

        props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

        qpas (list): List of {q,p,a,s} MATRIX objects for surface 1.

        coeff (CMATRIX(ntraj x 1)): The complex coefficient matrix for the TBF on all surfaces.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

    Returns:
        qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

        qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

        b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

        b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).

    """

    dt = params["dt"]
    decpl = params['decpl_den'] # 
    beta = params['linfit_beta']
    iM = params["iM"]  # MATRIX(ndof, 1)


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
    nstates = len(set(surf_ids))

    x_dofs = list(range(ndof))
    invM = MATRIX(ndof, ndof)
    for i in x_dofs:
        invM.set(i, i,  iM.get(i, 0))

    # The new set (propagated)
    q_new = MATRIX(q_old)
    p_new = MATRIX(p_old)
    a_new = MATRIX(a_old)
    s_new = MATRIX(s_old)

    unsorted_pairs = []
    for i in range(nstates):
        unsorted_pairs.append([i, surf_pops[i]])

    sorted_pairs = merge_sort(unsorted_pairs) 

    sorted_states = [0]*nstates
    for i in range(nstates):
        sorted_states[nstates-1-i] = sorted_pairs[i][0]

    #sorted_pops = sorted(surf_pops, reverse = True)
    #sorted_states = []
    #for n in range(nstates):
    #    sorted_states.append(surf_pops.index(sorted_pops[n]))

   
    # The properties for the trajectories on the most populated surface - reference for synchronizing
    q_new_on_surf_ref = None
    p_new_on_surf_ref = None
    a_new_on_surf_ref = None
    s_new_on_surf_ref = None

#    iref_pop_state = True
    for nindex, n in enumerate(sorted_states): # over all states, but starting with the most populated one

        traj_on_surf = [index for index, traj_id in enumerate(surf_ids) if traj_id == n]
        ntraj_on_surf = len(traj_on_surf)  

        q_on_surf = MATRIX(ndof, ntraj_on_surf)
        p_on_surf = MATRIX(ndof, ntraj_on_surf)
        a_on_surf = MATRIX(ndof, ntraj_on_surf)
        s_on_surf = MATRIX(ndof, ntraj_on_surf)

        pop_submatrix(q_old, q_on_surf, x_dofs, traj_on_surf)
        pop_submatrix(p_old, p_on_surf, x_dofs, traj_on_surf)
        pop_submatrix(a_old, a_on_surf, x_dofs, traj_on_surf)
        pop_submatrix(s_old, s_on_surf, x_dofs, traj_on_surf)


        # "Independent" evolution of trajectories
        if surf_pops[n] > decpl:
            coeff_on_surf = CMATRIX(ntraj_on_surf,1)  # coefficients for the TBFs on the surface n
            pop_submatrix(coeff, coeff_on_surf, traj_on_surf, [0])

            qpas_on_surf = [q_on_surf, p_on_surf, a_on_surf, s_on_surf] # variables for the TBFs on the surface n
            mom, r, gmom, gr = qtag_mom.momentum(params, ndof, ntraj_on_surf, qpas_on_surf, coeff_on_surf)

            # mom - MATRIX(ndof, ntraj_on_surf) - quantum momentum for q -  Im( nubla_{\alp} \psi / \psi)
            # r - MATRIX(ndof, ntraj_on_surf) - quantum momentum for s   -  Re( nubla_{\alp} \psi / \psi)
            # gmom - MATRIX(ndof, ntraj_on_surf) - gradient of the fit of mom ~ nubla_{\alp} Im( nubla_{\alp} \psi / \psi) - to update alphas
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
            if q_sync_method == 1:
                push_submatrix(q_new, q_new_on_surf_ref, x_dofs, traj_on_surf)
            if p_sync_method == 1:
                push_submatrix(p_new, p_new_on_surf_ref, x_dofs, traj_on_surf)
            if a_sync_method == 1:
                push_submatrix(a_new, a_new_on_surf_ref, x_dofs, traj_on_surf)
            if s_sync_method == 1:
                push_submatrix(s_new, s_new_on_surf_ref, x_dofs, traj_on_surf)


    qpas_new = [q_new, p_new, a_new, s_new, surf_ids]

    ov_no = qtag_calc.new_old_overlap(ndof, ntraj, nstates, qpas, qpas_new)
    btot = ov_no*coeff

    return qpas_new, btot



def cls_force(univ,mss,mom_calc,props,model_params,qpas1,c1_new,qpas2,c2_new,norm2,beta):
    """Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection vectors *b1* and *b2*, where the motion of both sets of functions are calculated via classical forces computed at their centers. Also necessary are the functions for calculating momentum (*mom_calc*) and basis updates (*props*).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        mss (dictionary): Dictionary containing multi-surface scheme parameters.

        mom_calc (function object): The function object needed to calculate the momentum, as defined by qtag_config.

        props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

        model_params (dictionary): Dictionary containing the potential parameters.

        qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

        c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

        c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

        norm2 (float): The population on surface 2.

        beta (float): Parameter determining tolerance in the momentum linear fitting algorithm.

    Returns:
        qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

        qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

        b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

        b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
    """

    ndof,ntraj=univ['ndof'],univ['ntraj']
    decpl=mss['decpl']
    qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3];cls_prop=props[4]
    qvals,pvals,avals,svals=MATRIX(qpas1[0]),MATRIX(qpas1[1]),MATRIX(qpas1[2]),MATRIX(qpas1[3])
    qvalsn,pvalsn,avalsn,svalsn=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)

    mom,r,gmom,gr=mom_calc(univ,beta,qpas1,c1_new)

    for i in range(ndof):
        qn=qprop(univ,i,qvals.row(i),mom.row(i));pn=pprop(mom.row(i));an=aprop(univ,i,avals.row(i),gmom.row(i));sn=sprop(univ,i,svals.row(i))
        for j in range(ntraj):
            qvalsn.set(i,j,qn.get(j))
            pvalsn.set(i,j,pn.get(j))
            avalsn.set(i,j,an.get(j))
            svalsn.set(i,j,sn.get(j))

    qpas1n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]
    ov_no=qtag_calc.overlap(ntraj,qpas1n,qpas1)
    b1=ov_no*c1_new

    if norm2 < decpl:
        qvals,pvals=MATRIX(qpas2[0]),MATRIX(qpas2[1])
        qvals_cls,pvals_cls=cls_prop(univ,qvals,pvals,model_params)
        qpas2n=[MATRIX(qvals_cls),MATRIX(pvals_cls),MATRIX(avalsn),MATRIX(svalsn)]
    else:
        qvals,pvals,avals,svals=MATRIX(qpas2[0]),MATRIX(qpas2[1]),MATRIX(qpas2[2]),MATRIX(qpas2[3])
        mom,r,gmom,gr=mom_calc(univ,beta,qpas2,c2_new)
        for i in range(ndof):
            qn2=qprop(univ,i,qvals.row(i),mom.row(i));pn2=pprop(mom.row(i));an2=aprop(univ,i,avals.row(i),gmom.row(i));sn2=sprop(univ,i,svals.row(i))
            for j in range(ntraj):
                qvalsn.set(i,j,qn2.get(i,j))
                pvalsn.set(i,j,pn2.get(i,j))
                avalsn.set(i,j,an2.get(i,j))
                svalsn.set(i,j,sn2.get(i,j))
        qpas2n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]

    ov_no=qtag_calc.overlap(ntraj,qpas2n,qpas2)
    b2=ov_no*c2_new

    return(qpas1n,qpas2n,b1,b2)

def fixed():
    return ()

def mean_field(univ,mss,mom_calc,props,model_params,qpas1,c1_new,qpas2,c2_new,norm2,beta):
    """Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), 
       as well as their corresponding projection vectors *b1* and *b2*, where the basis parameters are calculated 
       according to an "average" surface, defined by the average momentum calculed in the *mom_avg* function. 
       Also necessary are the functions for calculating momentum (*mom_calc*, although specifically mom_avg here) 
       and basis updates (*props*).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        mss (dictionary): Dictionary containing multi-surface scheme parameters.

        mom_calc (function object): The function object needed to calculate the momentum, as defined by qtag_config.

        props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

        model_params (dictionary): Dictionary containing potential parameters.

        qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

        c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

        c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

        norm2 (float): The population on surface 2

        beta (float): Parameter determining tolerance in the momentum linear fitting algorithm.

    Returns:
        qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

        qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

        b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

        b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
    """

    ndof,ntraj=univ['ndof'],univ['ntraj']
    decpl=mss['decpl']
    qvalsn,pvalsn,avalsn,svalsn=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)
    qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]
    qvals,pvals,avals,svals=MATRIX(qpas1[0]),MATRIX(qpas1[1]),MATRIX(qpas1[2]),MATRIX(qpas1[3])

    mom,r,gmom,gr=mom_calc(univ, beta, qpas1, c1_new, qpas2, c2_new)
    for i in range(ndof):
        qn=qprop(univ,i,qvals.row(i),mom.row(i),mss);pn=pprop(mom.row(i));an=aprop(univ,i,avals.row(i),gmom.row(i),mss);sn=sprop(univ,i,svals.row(i))
        for j in range(ntraj):
            qvalsn.set(i,j,qn.get(i,j))
            pvalsn.set(i,j,pn.get(i,j))
            avalsn.set(i,j,an.get(i,j))
            svalsn.set(i,j,sn.get(i,j))

    qpas1n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]
    qpas2n=[MATRIX(qvalsn),MATRIX(pvalsn),MATRIX(avalsn),MATRIX(svalsn)]

    ov_no=qtag_calc.overlap(ntraj,qpas1n,qpas1)
    b1=ov_no*c1_new

    ov_no=qtag_calc.overlap(ntraj,qpas2n,qpas2)
    b2=ov_no*c2_new

    return(qpas1n,qpas2n,b1,b2)


def two_surf(univ,mss,mom_calc,props,model_params,qpas1,c1_new,qpas2,c2_new,norm2,beta):
    """Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), 
       as well as their corresponding projection vectors *b1* and *b2*, where the basis parameters are 
       calculated in pairs, although the surface assignment alternates with each function (i.e. even basis 
       are ground surface, odd basis are excited surface). Also necessary are the functions for calculating 
       momentum (*mom_calc*, although specifically mom_avg here) and basis updates (*props*).

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        mss (dictionary): Dictionary containing multi-surface scheme parameters.

        mom_calc (function object): The function object needed to calculate the momentum, as defined by qtag_config.

        props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

        model_params (dictionary): Dictionary containing potential parameters.

        qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

        c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

        c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

        norm2 (float): The population on surface 2

        beta (float): Parameter determining tolerance in the momentum linear fitting algorithm.

    Returns:
        qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

        qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

        b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

        b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).

    """

    ndof,ntraj=univ['ndof'],univ['ntraj']
    qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]

    qvals1,pvals1,avals1,svals1=MATRIX(qpas1[0]),MATRIX(qpas1[1]),MATRIX(qpas1[2]),MATRIX(qpas1[3])
    qvals2,pvals2,avals2,svals2=MATRIX(qpas2[0]),MATRIX(qpas2[1]),MATRIX(qpas2[2]),MATRIX(qpas2[3])
    qvals1n,pvals1n,avals1n,svals1n=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)
    qvals2n,pvals2n,avals2n,svals2n=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)
    b1,b2=CMATRIX(ntraj,1),CMATRIX(ntraj,1)
    dummy=CMATRIX(ntraj,ntraj)

    mom1,r1,gmom1,gr1=mom_calc(univ,beta,qpas1,c1_new)
    mom2,r2,gmom2,gr2=mom_calc(univ,beta,qpas2,c2_new)

    for i in range(ndof):
        qn1=qprop(univ,i,qvals1.row(i),mom1.row(i),mss);pn1=pprop(mom1.row(i));an1=aprop(univ,i,avals1.row(i),gmom1.row(i),mss);sn1=sprop(univ,i,svals1.row(i))
        qn2=qprop(univ,i,qvals2.row(i),mom2.row(i),mss);pn2=pprop(mom2.row(i));an2=aprop(univ,i,avals2.row(i),gmom2.row(i),mss);sn2=sprop(univ,i,svals2.row(i))
        for j in range(ntraj):
            if (j%2==0):
                qvals1n.set(i,j,qn1.get(j)); qvals2n.set(i,j,qn1.get(j))
                pvals1n.set(i,j,pn1.get(j)); pvals2n.set(i,j,pn1.get(j))
                avals1n.set(i,j,an1.get(j)); avals2n.set(i,j,an1.get(j))
                svals1n.set(i,j,sn1.get(j)); svals2n.set(i,j,sn1.get(j))
            else:
                qvals1n.set(i,j,qn2.get(j)); qvals2n.set(i,j,qn2.get(j))
                pvals1n.set(i,j,pn2.get(j)); pvals2n.set(i,j,pn2.get(j))
                avals1n.set(i,j,an2.get(j)); avals2n.set(i,j,an2.get(j))
                svals1n.set(i,j,sn2.get(j)); svals2n.set(i,j,sn2.get(j))

    qpas1n=[MATRIX(qvals1n),MATRIX(pvals1n),MATRIX(avals1n),MATRIX(svals1n)]
    qpas2n=[MATRIX(qvals2n),MATRIX(pvals2n),MATRIX(avals2n),MATRIX(svals2n)]

    ov_no=qtag_calc.overlap(ntraj,qpas1n,qpas1)
    b1=ov_no*c1_new

    ov_no=qtag_calc.overlap(ntraj,qpas2n,qpas2)
    b2=ov_no*c2_new
	
    return(qpas1n,qpas2n,b1,b2)
