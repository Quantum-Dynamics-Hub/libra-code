#*********************************************************************************                     
#* Copyright (C) 2022-2023 Matthew Dutra and Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: compute
   :platform: Unix, Windows
   :synopsis: This module implements functions for doing QTAG dynamics
       List of functions:  
           * time_overlap(nQ, nP, nA, nS, nstate, oQ, oP, oA, oS, ostate)
           * qtag_momentum(dyn_params,qpas,coeff_on_surf)
           * propagate_basis(_q, _p, _alp, _s, _states, coeff, dyn_params, surf_pops)
           * qtag_pops(surf_ids, coeff, S, target_states)
           * qtag_energy(coeff, H)
           * psi(ndof, ntraj_on_surf, qpas, c, x0)
           * wf_calc_nD(dyn_params, plt_params, prefix)
           * run_qtag(_q, _p, _alp, _s, _states, _coeff, _iM, _dyn_params, _compute_model, _model_params)

.. moduleauthor:: Matthew Dutra and Alexey V. Akimov

"""

__author__ = "Matthew Dutra, Alexey V. Akimov"
__copyright__ = "Copyright 2022 Matthew Dutra, Alexey V. Akimov"
__credits__ = ["Matthew Dutra, Alexey V. Akimov"]
__license__ = "GNU-3"
__version__ = "1.0"
__maintainer__ = "Alexey V. Akimov"
__email__ = "alexvakimov@gmail.com"
__url__ = "https://quantum-dynamics-hub.github.io/libra/index.html"


import os
import sys
import math
import copy
import time
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import util.libutil as comn
import libra_py.units as units
import libra_py.data_outs as data_outs
import libra_py.data_savers as data_savers

from . import save


def time_overlap(nQ, nP, nA, nS, nstate, oQ, oP, oA, oS, ostate):
    """
    Computes the time-overlap <G_new|G_old>

    n - new 
    o - old
    
    Args: 
    
        nQ (MATRIX(ndof, ntraj)): Matrix containing basis positions after a timestep dt
        
        nP (MATRIX(ndof, ntraj)): Matrix containing basis momenta after a timestep dt
        
        nA (MATRIX(ndof, ntraj)): Matrix containing basis widths after a timestep dt
        
        nS (MATRIX(ndof, ntraj)): Matrix containing basis phases after a timestep dt
        
        nstate (list of ints): The list specifying which electronic state each updated GBF belongs to
        
        oQ (MATRIX(ndof, ntraj)): Matrix containing basis positions before a timestep dt
        
        oP (MATRIX(ndof, ntraj)): Matrix containing basis momenta before a timestep dt
        
        oA (MATRIX(ndof, ntraj)): Matrix containing basis widths before a timestep dt
        
        oS (MATRIX(ndof, ntraj)): Matrix containing basis phases before a timestep dt
        
        ostate (list of ints): The list specifying which electronic state each old GBF belongs to
        
    Returns:
    
        st (CMATRIX(ntraj, ntraj)): The time-overlap matrix between sets of GBFs across a timestep dt
    """

    ndof = nQ.num_of_rows
    ntraj = nQ.num_of_cols
    nstates = len(set(nstate))

    nA_half = 0.5*nA
    oA_half = 0.5*oA    

    St = gwp_overlap_matrix(nQ, nP, nS, nA_half, Py2Cpp_int(nstate), oQ, oP, oS, oA_half, Py2Cpp_int(ostate))

    return St


def qtag_momentum(_q, _p, _a, _s, coeff_on_surf, dyn_params):
    """Calculates the single-surface momentum for a set of basis functions, received from
       the `propagate( )' function.

    Args:

        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`d_weight`]** (int) : parameter used to specify whether the fitted
              momentum values should be density-weighted in the linear fitting function. It
              is highly recommended to leave this value as 1 (on).  [ default: 1 ]

          * **dyn_params[`linfit_beta`]** (float) : parameter used to specify the convergence
              criterion for the linear fitting algorithm. A smaller value indicates a
              stricter fit, although convergence issues may arise when going below 1e-5  
              [ default: 1e-1 ]

        _q ( MATRIX(ndof, ntraj_on_surf) ): coordinates of the "classical" particles [units: Bohr]
        
        _p ( MATRIX(ndof, ntraj_on_surf) ): momenta of the "classical" particles [units: a.u. of momenta]
        
        _a ( MATRIX(ndof, ntraj_on_surf) ): widths of the GBFs [ units: Bohr^-1 ]
        
        _s ( MATRIX(ndof, ntraj_on_surf) ): phases of all GBFs [ units: no ]

        coeff_on_surf (CMATRIX(ntraj_on_surf x 1)): The complex coefficient matrix for the TBF
        on the relevant surface.

    Returns:
        mom (MATRIX(ndof, ntraj)): Matrix containing trajectory momenta.

        r (MATRIX(ndof, ntraj)): Matrix containing the real complement to the trajectory momenta.

        gmom (MATRIX(ndof, ntraj)): Matrix containing the momentum gradients at each trajectory
        location.

        gr (MATRIX(ndof, ntraj)): Matrix containing the real complement to *gmom*.
    """

    params = dict(dyn_params)

    critical_params = [ ]
    default_params = { "d_weight":1, "linfit_beta":1e-1, "mom_calc_type":1 }
    comn.check_input(params, default_params, critical_params)
 
    #Extract parameters from params dict...
    mom_calc_type = params['mom_calc_type']
    beta = params['linfit_beta']
    d_weight = params['d_weight']

    ndof = _q.num_of_rows
    ntraj_on_surf = _q.num_of_cols

    #Assign MATRIX and CMATRIX objects...
    mom = MATRIX(ndof, ntraj_on_surf)
    r = MATRIX(ndof, ntraj_on_surf)
    gmom = MATRIX(ndof, ntraj_on_surf)
    gr = MATRIX(ndof, ntraj_on_surf)

    psi_tot = CMATRIX(ntraj_on_surf,1)
    grad_psi = CMATRIX(ndof,1)
    deriv_term = CMATRIX(ndof,1)

    q_on_surf = MATRIX(_q)
    p_on_surf = MATRIX(_p)
    a_on_surf = MATRIX(_a)
    s_on_surf = MATRIX(_s)

    #Calculate momentum for each trajectory (i) as Im(grad(psi)/psi)...
    for i in range(ntraj_on_surf):
        psi_sum = complex(0.0,0.0)
        for dof in range(ndof):
            grad_psi.set(dof,0,0+0j)

        for j in range(ntraj_on_surf):
            pre_exp_full = 1.0
            exp_full = 1.0
            c2 = coeff_on_surf.get(j)

            for dof in range(ndof):
                q1 = q_on_surf.get(dof,i)
                q2 = q_on_surf.get(dof,j)
                dq = q1-q2

                p1 = p_on_surf.get(dof,i)
                p2 = p_on_surf.get(dof,j)

                a1 = a_on_surf.get(dof,i)
                a2 = a_on_surf.get(dof,j)

                s1 = s_on_surf.get(dof,i)
                s2 = s_on_surf.get(dof,j)

                pre_exp_full *= (a2/np.pi)**0.25
                exp_full *= np.exp(-0.5*a2*dq**2+1.0j*(p2*dq+s2))
                deriv_term.set(dof, complex(-a2*dq,p2) )

            psi_at_j = c2*pre_exp_full*exp_full
            psi_sum += psi_at_j

            for dof in range(ndof):
                grad_psi.set(dof, grad_psi.get(dof)+psi_at_j*deriv_term.get(dof))

        psi_tot.set(i, psi_sum)
        for dof in range(ndof):
            mom.set(dof,i,(grad_psi.get(dof)/psi_sum).imag)
            r.set(dof,i,(grad_psi.get(dof)/psi_sum).real)
            gmom.set(dof,i,0.0)
            gr.set(dof,i,0.0)

    #For linear fitting of momentum, procedure is least squares fitting of type Ax=B...
    if mom_calc_type == 1:
        for dof in range(ndof):
            A=MATRIX(2,2);B=MATRIX(2,2);x=MATRIX(2,2)

            for m in range(2):
                for n in range(2):
                    elem_A = 0; elem_B = 0

                    for i in range(ntraj_on_surf):
                        q = q_on_surf.get(dof,i)
                        pimag = mom.get(dof,i)
                        preal = r.get(dof,i)

                        if d_weight == 1:
                            z = psi_tot.get(i)
                            zstar = np.conj(z)
                        else:
                            z=1+0j; zstar=1-0j

                        elem_A += q**(m+n)*(z*zstar).real
                        if n == 0:
                            elem_B += pimag*q**(m)*(z*zstar).real
                        elif n == 1:
                            elem_B += preal*q**(m)*(z*zstar).real

                    A.set(m,n,elem_A)
                    B.set(m,n,elem_B)

            solve_linsys(A,B,x,beta,200000)
            for i in range(ntraj_on_surf):
                q = q_on_surf.get(dof,i)
                aa=x.get(0,0)+x.get(1,0)*q
                bb=x.get(0,1)+x.get(1,1)*q

                mom.set(dof,i,aa); r.set(dof,i,bb)
                gmom.set(dof,i,x.get(1,0)); gr.set(dof,i,x.get(1,1))

    return mom, r, gmom, gr



def propagate_basis(_q, _p, _alp, _s, _states, coeff, dyn_params, surf_pops):
    """Computes the updated basis parameters q, p, a, and s according to a specified single-surface and multi-surface
    calculation scheme, as defined by the update_method and sync_method sets of parameters in dyn_params.

    Args:

        dyn_params (dict): Dictionary containing simulation parameters.
 
          * **dyn_params[`decpl_den`]** (float) : a parameter controlling independence of trajectories. If the population
              on a given state is larger that this parameter, the trajectories evolve according to their own
              quantum momentum. If the population is less than this number, the trajectories evolve according
              to the quantum momentum for the most populated state.  [ default: 0.25 ]

        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        
        _p ( MATRIX(nnucl, ntraj) ): momenta of the "classical" particles [units: a.u. of momenta]
        
        _alp ( MATRIX(nnucl, ntraj) ): widths of the GBFs [ units: Bohr^-1 ]
        
        _s ( MATRIX(1, ntraj) ): phases of all GBFs [ units: no ]
        
        _states ( intList, or list of ntraj ints ): the quantum state of each trajectory

        coeff (CMATRIX(ntraj, 1)): The complex coefficient matrix for the TBF on all surfaces.

        surf_pops (list of floats): List of surface populations.

    Returns:
    
        q_new (MATRIX(ndof, ntraj)): the updated matrix object containing the basis positions
        
        p_new (MATRIX(ndof, ntraj)): the updated matrix object containing the basis momenta
        
        a_new (MATRIX(ndof, ntraj)): the updated matrix object containing the basis widths
        
        s_new (MATRIX(ndof, ntraj)): the updated matrix object containing the basis phases 
        
        surf_ids (list of ints): the list specifying which electronic state each trajectory belongs to
        
        coeff (CMATRIX(ntraj, 1)):  the updated complex coefficient matrix for the TBF on all surfaces
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
    q_old = MATRIX(_q)
    p_old = MATRIX(_p)
    a_old = MATRIX(_alp)
    s_old = MATRIX(_s)    
    surf_ids = _states

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

#            qpas_on_surf = [q_on_surf, p_on_surf, a_on_surf, s_on_surf] # variables for the TBFs on the surface n
            mom, r, gmom, gr = qtag_momentum(q_on_surf, p_on_surf, a_on_surf, s_on_surf, coeff_on_surf, dyn_params)

            # mom - MATRIX(ndof, ntraj_on_surf) - quantum momentum for q -  Im( nabla_{\alp} \psi / \psi)
            # r - MATRIX(ndof, ntraj_on_surf) - quantum momentum for s   -  Re( nabla_{\alp} \psi / \psi)
            # gmom - MATRIX(ndof, ntraj_on_surf) - gradient of the fit of mom ~ nabla_{\alp} Im( nubla_{\alp} \psi / \psi) - to update alphas
            # gr - MATRIX(ndof, ntraj_on_surf) - gradient of the fit of r ~ nubla_{\alp} Re( nubla_{\alp} \psi / \psi)  - to update alphas

            # mom_tmp = momentum(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff);
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


    #========= Re-expansion of coeffients =====================
    a_new_half = 0.5* a_new    
    GG = gwp_overlap_matrix(q_new, p_new, s_new, a_new_half, Py2Cpp_int(surf_ids), 
                            q_new, p_new, s_new, a_new_half, Py2Cpp_int(surf_ids))

    tmp = FullPivLU_rank_invertible(GG)
    if tmp[1]==0:
        print("GBF overlap matrix is not invertible.\n Exiting...\n")
        sys.exit(0)

    invGG = CMATRIX(ntraj, ntraj)
    FullPivLU_inverse(GG, invGG)


    st = time_overlap(q_new, p_new, a_new, s_new, surf_ids, q_old, p_old, a_old, s_old, surf_ids )
    coeff = st * coeff
    coeff = invGG * coeff

    return q_new, p_new, a_new, s_new, surf_ids , coeff




def qtag_pops(surf_ids, coeff, S, target_states):
    """Returns a single-surface population *n*, calculated from the single-surface basis coefficients *c* 
       and the appropriate overlap matrix *ov* as <c^T|ov|c>. Note that this is equivalent to the norm for 
       a single-surface system.

    Args:
        surf_ids (list of ints): List containing the trajectory indices on various states.

        coeff ( CMATRIX(ntraj, 1) ): The GBF expansion coefficients (in the super-basis) 

        S (CMATRIX(ntraj, ntraj) ): Basis functions super-overlap

        target_states (list of ints): List of states for which the norm should be calculated.

    Returns:
        pops (list of floats): Surface population. The imaginary part should be zero.

    """

    pops = []
    for n in target_states:
        # Extract indices of the trajectories that sit on a given state n
        traj_on_surf = [index for index, traj_id in enumerate(surf_ids) if traj_id == n]
        ntraj_on_surf = len(traj_on_surf)

        # Extract the overlap matrix of the GBFs on that surface
        ov_surf = CMATRIX(ntraj_on_surf,ntraj_on_surf)
        pop_submatrix(S,     ov_surf, traj_on_surf, traj_on_surf)

        # Extract the coefficients of the GBFs on that surface
        c_surf = CMATRIX(ntraj_on_surf,1)
        pop_submatrix(coeff, c_surf,  traj_on_surf, [0])

        pops.append((c_surf.H()*ov_surf*c_surf).get(0).real)

    return pops


def qtag_energy(coeff, H):
    """Returns the system energy *e*, calculated from the total basis coefficients *coeff* 
       and the total Hamiltonian *H* as <G|H|G>.
       E = C^+ * H * C

    Args:
        coeff (CMATRIX(ntraj, 1) ): The ntraj-by-1 complex matrix of basis coefficients.

        H (CMATRIX(ntraj, ntraj) ): The full system super-Hamiltonian (all trajectories and surfaces).

    Returns:
        e (float): Total energy of the system. The imaginary part should be zero.
    """

    e = (coeff.H() * H * coeff).get(0).real

    return e


def wf_calc_nD(_params):
    """Computes the wavefunction on an nD-dimensional grid, as specified by the parameters obtained
    from the dyn_params and plt_params dictionaries.

    Args:
        _params (dict): Dictionary containing control parameters.

          * **_params[`states`]** (list of ints) : list of all states in data

          * **_params[`grid_dims`]** (list of floats) : the total number of basis functions to be
              placed on each surface. For Gaussian, the list has only one element, specifying the
              number of basis functions per surface. Note that the total number of basis functions
              will then be *nstates*-by-*prod(grid_dims)*

          * **_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]
                  
          * **_params[`xmin`]** (list of floats) : list of minima in the independent coordinate space
              at which the wavefunction should be calculated
              
          * **_params[`xmax`]** (list of floats) : list of maxima in the independent coordinate space
              at which the wavefunction should be calculated
              
          * **_params[`npoints`]** (list of ints) : the number of points used to compute the wavefunction
              along each coordinate (i.e. the calculation grid)
              
          * **_params[`snaps`]** (list of ints) : the snapshots at which to calculate the wavefunction

        prefix (str): The name of the directory containing the q, p, a, and s trajectory data.
    """


    params = dict(_params)

    critical_params = [ "grid_dims" ]
    default_params = { "states":[0], "ndof":1, "xmin":[-4.0], "xmax":[4.0],
                       "npoints":[100], "snaps":[0], "prefix":"out"
                     }
    comn.check_input(params, default_params, critical_params)


    #Collect simulation parameters from dyn_params dict...
    #nsteps = params["nsteps"]
    states = params["states"]
    grid_dims = params["grid_dims"]
    ndof = params["ndof"]
    prefix = params["prefix"]

    nstates = len(states)
    xmin = params["xmin"]
    xmax = params["xmax"]
    npoints = params["npoints"]

    #This is standard Libra convention for wf files...
    data_type1="wfcr"; data_type2="dens"; data_type3="rep_0"

    #Determine ntraj from dyn_params grid_dims variable...
    ntraj = 1   # this is the number of trajectories per surface
    for i in range(len(grid_dims)):
        ntraj *= grid_dims[i]
        #print(F"i = {i},  grid_dims[{i}] = {grid_dims[i]}, ntraj = {ntraj}")

    #Define which timesteps to calculate the wf for...
    snaps = params['snaps']

    #Open directory with coeffs, q, p, a, s data...
    qfile = open(prefix+"/q.txt")
    pfile = open(prefix+"/p.txt")
    afile = open(prefix+"/a.txt")
    cfile = open(prefix+"/coeffs.txt")

    #Make the output directory...
    if not os.path.isdir(prefix+"/wfc"):
        os.mkdir(prefix+"/wfc")

    #Create the mesh to calculate the wf on...
    grid_bounds = np.mgrid[tuple(slice(xmin[dof],xmax[dof],complex(0,npoints[dof])) for dof in range(ndof))]

    # How many grid points in total. npoints[0], npoints[1], ... are the numbers of the points
    # along dimensions 0, 1, etc.
    elems = 1
    for i in npoints:
        elems *= i

    b = grid_bounds.flatten()

    # Coordinates of all the grid points
    wfpts = []
    for i in range(elems):
        index = i
        coords = []

        while index < len(b):
            coords.append(b[index])
            index += elems
        wfpts.append(coords)

    #Read the trajectory data and compute the wf on the mesh...
    Qdata = qfile.readlines()
    Pdata = pfile.readlines()
    Adata = afile.readlines()
    Coeff = cfile.readlines()

    for isnap in snaps:
        qdata = Qdata[isnap].strip().split()
        pdata = Pdata[isnap].strip().split()
        adata = Adata[isnap].strip().split()
        coeffs = Coeff[isnap].strip().split()

        #print(len(coeffs), F"nstates = {nstates}, ntraj = {ntraj}  isnap = {isnap}")
        outfile = open(F"{prefix}/wfc/{data_type1}_snap_{isnap}_{data_type2}_{data_type3}","w")

        for pt in wfpts:
            for nn in pt:
                outfile.write(str(nn)+" ")

            idata = 0
            for state in range(nstates):
                wf = 0+0j
                for j in range(ntraj):
                    re_indx = 2 * ntraj * state + 2 * j + 0
                    im_indx = 2 * ntraj * state + 2 * j + 1

                    #print(F" state = {state}, traj = {j}, re_indx = {re_indx}, im_indx = {im_indx}")

                    re = float(coeffs[re_indx])
                    im = float(coeffs[im_indx])

                    coeff = complex(re, im)
             
                    #coeff = complex(float(coeffs[2*(state*ntraj+j) ]), float(coeffs[2*(state*ntraj+j)+1])) 

                    gaus = 1.0+0.0j
                    for dof in range(ndof):
                        qj = float(qdata[idata])
                        pj = float(pdata[idata])
                        aj = float(adata[idata])
                        idata +=1

                        q = pt[dof]
                        gaus*=(aj/np.pi)**0.25 * np.exp(-0.5*aj*(q-qj)**2 + 1j*pj*(q-qj)  )

                    wf += coeff*gaus
                outfile.write(str(abs(wf)**2)+" ")
            outfile.write("\n")


def run_qtag(_q, _p, _alp, _s, _states, _coeff, _iM, _dyn_params, _compute_model, _model_params):
    """
    Args: 
        _q ( MATRIX(nnucl, ntraj) ): coordinates of the "classical" particles [units: Bohr]
        _p ( MATRIX(nnucl, ntraj) ): momenta of the "classical" particles [units: a.u. of momenta]
        _alp ( MATRIX(nnucl, ntraj) ): widths of the GBFs [ units: Bohr^-1 ]
        _s ( MATRIX(1, ntraj) ): phases of all GBFs [ units: no ]
        _states ( intList, or list of ntraj ints ): the quantum state of each trajectory
        _coeff ( CMATRIX(nstates x ntraj) ): amplitudes of the GBFs
        _iM ( MATRIX(nnucl, 1) ): masses of classical particles [units: a.u.^-1]
        _dyn_params ( dictionary ): parameters controlling the execution of the dynamics
            Can contain:

            * **q_update_method** (int): the parameter specifying whether positions should be updated (adaptable) or frozen  

                - 0: frozen
                - 1: adaptable [ default ]


            * **p_update_method** (int): the parameter specifying whether momenta should be updated (adaptable) or frozen  

                - 0: frozen
                - 1: adaptable [ default ]


            * **a_update_method** (int): the parameter specifying whether widths should be updated (adaptable) or frozen  

                - 0: frozen [ default ]
                - 1: adaptable


            * **s_update_method** (int): the parameter specifying whether phases should be updated (adaptable) or frozen 

                - 0: frozen [ default ]
                - 1: adaptable (not available yet)


            * **q_sync_method** (int): the parameter specifying whether positions should be synchronized on low pop surfaces  

               - 0: unsynchronized (leave at initial value)
               - 1: synchronized to most populated state [ default ]


            * **p_sync_method** (int): the parameter specifying whether momenta should be synchronized on low pop surfaces  

              - 0: unsynchronized (leave at initial value)
              - 1: synchronized to most populated state [ default ]


            * **a_sync_method** (int): the parameter specifying whether widths should be synchronized on low pop surfaces  

              - 0: unsynchronized (leave at initial value) [ default ]
              - 1: synchronized to most populated state 


            * **s_sync_method** (int): the parameter specifying whether phases should be synchronized on low pop surfaces  

              - 0: unsynchronized (leave at initial value) [ default ]


            * **decpl_den** (float): a critical population above which a surface can evolve trajectories independently  


            * **qtag_pot_approx_method** (int): the method for approximating potential matrix elements in a basis  

              - 0: BAT (Bra-ket Averaged Taylor expansion)
              - 1: LHA (Local Harmonic Approximation)
              - 2: LHAe (Local Harmonic Approximation w/ exact Gaussian coupling elements)
              - 3: BATe (Bra-ket Averaged Taylor expansion w/ exact Gaussian coupling elements)


            * **mom_calc_type** (int): how to modify the raw basis momenta for numerical stability  

              - 0: unmodified (note that these are frequently unstable!)
              - 1: linear fitting of single-surface momenta


            * **linfit_beta** (float): a threshold for the linear fitting algorithm when needed for momenta modification  

 
            * **prefix** (string): name of the directory where all the results will be stored [ default: "out" ]


            * **hdf5_output_level** (int): the level of output for HDF5 printing of data [ default: -1 ] 


            * **txt2_output_level** (int): the level of output for txt printing of data; 3 is the current maximum  [ default: 0 ]


            * **progress_frequency** (int): how often to print out the results [ default: 1 - every timestep] 


            * **properties_to_save** (list of strings): a list containing the desired output quantities 


            * **target_states** (list of ints): the indices of quantum states for which to compute populations


            * **dt** (float): the timestep between each iteration [ in a.u. of time ]


            * **nsteps** (int): the number of total iterations; this together with `dt` determines the simulated time 


        _compute_model ( PyObject ): the function that computes the Hamiltonian object
        _model_params ( dict ): parameters that are passed to the Hamiltonian-computing function

    """

    dyn_params = dict(_dyn_params)
    Q = MATRIX(_q)
    P = MATRIX(_p)
    A = MATRIX(_alp)
    S = MATRIX(_s)
    C = CMATRIX(_coeff)

    default_params = {
        "hdf5_output_level":-1, "prefix":"out", "use_compression":0, "compression_level":[0,0,0], 
        "mem_output_level":4, "txt2_output_level":0, "properties_to_save": [], "progress_frequency": 1,
        "target_states":[0]
    }
    critical_params = []
    comn.check_input(dyn_params, default_params, critical_params)

    ndof = Q.num_of_rows
    ntraj = Q.num_of_cols  # total number of trajectories 
    nstates = len( set(_states) )
    active_states = list(_states)

    #Rename variables locally for convenience...
    #states = sorted(_states)  # <---- TENTATIVELY comment
    #nstates = len(states)
    #active_state = dyn_params["active_state"]
    target_states = dyn_params["target_states"]
    nsteps = dyn_params["nsteps"]
    dt = dyn_params["dt"]
    iM = dyn_params["iM"]
    progress_frequency = dyn_params["progress_frequency"]
    prefix = dyn_params["prefix"]

    #Initialize savers...
    properties_to_save = dyn_params['properties_to_save']
    _savers = save.init_qtag_savers(dyn_params, _model_params, nsteps, ntraj, ndof, nstates)

    #Start simulation and walltime variables...
    walltime_start = time.time()
    t=0.0
    
    #Initialize the Hamiltonian objects in the C++ component of Libra...
    ham = nHamiltonian(nstates, nstates, ndof)
    ham.add_new_children(nstates, nstates, ndof, ntraj)    
    ham.init_all(_model_params["deriv_lvl"],1)
    _model_params.update({"timestep":0})
    
    ovlp = CMATRIX(ntraj, ntraj)
    hmat = CMATRIX(ntraj, ntraj)
    
    Coeff = CMATRIX(nstates,ntraj)

    # AVA: originally we use sorted states for some reason
    qtag_hamiltonian_and_overlap(Q, P, A, S, Coeff, Py2Cpp_int(active_states), iM, ham, 
                                 _compute_model, _model_params, dyn_params, ovlp, hmat)
    
    #Run the dynamics...
    for step in range(nsteps):
        # Built-in function for propagation in the non-othogonal basis
        propagate_electronic_qtag(0.5*dt, C, hmat, ovlp)
        
        #Calculate the total energy and surface populations...
        etot = qtag_energy(C, hmat)
        pops = qtag_pops(_states, C, ovlp, target_states)


        #Update the basis parameters according to the new wavefunction (ct_new)...
        Q, P, A, S, active_states, C = propagate_basis(Q, P, A, S, active_states, C, dyn_params, pops)

        #Compute the Hamiltonian and overlap matrix elements using the C++ routine 'qtag_ham_and_ovlp'...
        qtag_hamiltonian_and_overlap(Q, P, A, S, Coeff, Py2Cpp_int(active_states), iM, ham, 
                                     _compute_model, _model_params, dyn_params, ovlp, hmat)
    
        # Built-in function for propagation in the non-othogonal basis
        propagate_electronic_qtag(0.5*dt, C, hmat, ovlp)


        #Output the energy and populations to the notebook for the user to see...
        if step % progress_frequency == 0:
            print(etot, pops)
        
        #Save the specified data...
        save.save_qtag_data(_savers, dyn_params, step+1, etot, 0, pops, C, Q, P, A, S)        
        if _savers["txt2_saver"]!=None:
            _savers["txt2_saver"].save_data_txt( F"{prefix}", properties_to_save, "a", 0)


    #Print the total simulation time...
    walltime_end = time.time()
    print("Total wall time: ",walltime_end-walltime_start)

    return Q, P, A, S, active_states, C
    
