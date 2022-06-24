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
..module:: qtag_calc
  :platform: Unix, Windows
  :synopsis: This module contains "ground-level" functions for calculations, mostly output things like energy, norm, etc.

..moduleauthors :: Matthew Dutra, Alexey Akimov
"""

import os, sys
import numpy as np

from liblibra_core import *
from libra_py import data_outs

# CMATRIX qtag_psi(MATRIX& q, MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1, CMATRIX& Coeff);
def psi(ndof, ntraj_on_surf, qpas, c, x0):
    """Returns the (complex) wavefunction value *wf* at a given point *x0*, calculated using the single-surface basis parameters stored in *qpas* and coefficients *c*.

    Args:
        ndof (integer): Number of degrees of freedom.

        ntraj_on_surf (integer): Number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        c (CMATRIX): The ntraj_on_surf-by-1 complex matrix of basis coefficients.

        x0 (MATRIX): The matrix of coordinates [x_1,...,x_ndof] at which the wavefunction value should be calculated.

    Returns:
        wf (complex): Complex value of the wavefunction at x0.
    """

    qvals,pvals,avals,svals=qpas[0],qpas[1],qpas[2],qpas[3]

    wf=0+0j
    for i in range(ntraj_on_surf):
        prod=1.0
        for j in range(ndof):
            q1,p1,a1,s1=qvals.get(j,i),pvals.get(j,i),avals.get(j,i),svals.get(j,i)
            prod*=(a1/np.pi)**0.25*np.exp(-a1/2.0*(x0.get(j)-q1)**2+1j*(p1*(x0.get(j)-q1)+s1))
        wf+=c.get(i)*prod
    return(wf)

def energy(c,H):
    """Returns the system energy *e*, calculated from the total basis coefficients *c* and the total Hamiltonian *H* as <c^T|H|c>.

    Args:
        c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

        H (CMATRIX): The full system Hamiltonian (both surfaces + coupling).

    Returns:
        e (float): Energy. The imaginary part should be zero.
    """

    e=(c.H()*H*c).get(0).real
    return(e)

def norm(surf_ids,c,ov,states):
    """Returns a single-surface population *n*, calculated from the single-surface basis coefficients *c* and the appropriate overlap matrix *ov* as <c^T|ov|c>. Note that this is equivalent to the norm for a single-surface system.

    Args:
        surf_ids (list): List containing the trajectory indices on various states.

        c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

        ov (CMATRIX): The single-surface overlap matrix for which the population is to be calculated.

        states (list): List of states for which the norm should be calculated.

    Returns:
        pops (list of floats): Surface population. The imaginary part should be zero.
    """

    pops = []
    for n in states:
        traj_on_surf = [index for index, traj_id in enumerate(surf_ids) if traj_id == n]
        ntraj_on_surf = len(traj_on_surf)

        ov_surf = CMATRIX(ntraj_on_surf,ntraj_on_surf)
        c_surf = CMATRIX(ntraj_on_surf,1)

        pop_submatrix(ov,ov_surf,traj_on_surf,traj_on_surf)
        pop_submatrix(c,c_surf,traj_on_surf,[0])

        pops.append((c_surf.H()*ov_surf*c_surf).get(0).real)

    return(pops)


#def new_old_overlap(ndof,ntraj,states,qpaso,qpasn):
#    new_old_ov = CMATRIX(ntraj,ntraj)

def time_overlap(ndof,ntraj,states,qpasn,qpaso):
    """
    Computes the time-overlap <G_new|G_old>

    qpasn - new 
    qpaso - old
    """
    St = CMATRIX(ntraj,ntraj)

    for i in range(ntraj):
        for j in range(ntraj):
            new_old_ov.set(i,j,complex(0.0,0.0))

    #Extract the components of the qpas object into their constituent parts: q, p, a, s, surface IDs.
    qvals_old = MATRIX(qpaso[0]); qvals_new = MATRIX(qpasn[0])
    pvals_old = MATRIX(qpaso[1]); pvals_new = MATRIX(qpasn[1])
    avals_old = MATRIX(qpaso[2]); avals_new = MATRIX(qpasn[2])
    svals_old = MATRIX(qpaso[3]); svals_new = MATRIX(qpasn[3])
    surf_ids=qpasn[4]

    itot = 0; jtot = 0

    for n1 in states:

        #Extract the relevent trajectories for the n1-th state...
        traj_on_surf_n1 = [index for index, traj_id in enumerate(surf_ids) if traj_id == n1]
        ntraj_on_surf_n1 = len(traj_on_surf_n1)

        qvalso_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        pvalso_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        avalso_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        svalso_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)

        pop_submatrix(qvals_old,qvalso_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(pvals_old,pvalso_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(avals_old,avalso_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(svals_old,svalso_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)

        qvalsn_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        pvalsn_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        avalsn_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        svalsn_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)

        pop_submatrix(qvals_new,qvalsn_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(pvals_new,pvalsn_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(avals_new,avalsn_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(svals_new,svalsn_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)

        ii,jj = ntraj_on_surf_n1, ntraj_on_surf_n1
        n12_mat = CMATRIX(ii,jj)

        for i in range(ntraj_on_surf_n1):

            qi = qvalsn_surf_n1.col(i)
            pi = pvalsn_surf_n1.col(i)
            ai = avalsn_surf_n1.col(i)
            si = svalsn_surf_n1.col(i)

            for j in range(ntraj_on_surf_n1):

                qj = qvalso_surf_n1.col(j)
                pj = pvalso_surf_n1.col(j)
                aj = avalso_surf_n1.col(j)
                sj = svalso_surf_n1.col(j)

                n12_mat.set(i,j,gwp_overlap(qi,pi,si,ai/2,qj,pj,sj,aj/2))

        for i in range(ii):
            for j in range(jj):
                new_old_ov.set(i+itot,j+jtot,n12_mat.get(i,j))

        itot += ii
        jtot += jj

    return St


def basis_diag(m, dt, H, S, b):
    """Returns the updated basis coefficients for both surfaces, stored in a single vector *c_new*, 
       computed as `c_new=Z*exp(-i*dt*eps)*Z_dagger*b`. The variables eps and Z are the eigenvalues 
       and eigenvectors obtained from diagonalizing the full *m*-by-*m* Hamiltonian matrix *H* using the 
       solve_eigen internal function. Note that the projection vector *b* and the full overlap matrix *ov* are required.

    Args:
        m (integer): Dimension of the eigenvalue and eigenvector matrices. Usually equal to 2*ntraj.

        dt (real): The timestep parameter.

        H (CMATRIX): Total system Hamiltonian, dimension 2*ntraj-by-2*ntraj.

        ov (CMATRIX): Total overlap matrix, dimension 2*ntraj-by-2*ntraj.

        b (CMATRIX): Projection vector, dimensioned 2*ntraj-by-1.

    Returns:

        c_new (CMATRIX): Updated basis coefficient vector for both surfaces. Computed as c_new=Z*exp(-i*dt*eps)*Z_dagger*b.

    """

    evals = CMATRIX(m,m)
    evecs = CMATRIX(m,m)
    solve_eigen(H,ov,evals,evecs,0)
    ct = evecs.H()
    c_new = evecs*(exp_(evals,-dt*1.0j))*ct*b

    return(c_new)

def wf_calc_nD(dyn_params, plt_params, prefix):
    """Returns the initial basis parameters {q,p,a,s} as a list of  ndof-by-ntraj matrices *qpas*,
       based on the input contained in the dict *dyn_params*. The placement is randomly chosen from a
       Gaussian distribution centered about the wavepacket maximum with a standard deviation *rho_cut*,
       and each surface has the same number of trajectories. The corresponding Gaussian-distributed
       momenta are ordered so that the wavepacket spreads as x increases.

    Args:
        dyn_params (dict): Dictionary containing simulation parameters.

          * **dyn_params[`nsteps`]** (int) : the number of simulation steps

          * **dyn_params[`states`]** (int) : the list of states

          * **dyn_params[`grid_dims`]** (list of floats) : the total number of basis functions to be
              placed on each surface. For Gaussian, the list has only one element, specifying the
              number of basis functions per surface. Note that the total number of basis functions
              will then be *nstates*-by-*prod(grid_dims)*

          * **dyn_params[`ndof`]** (int) : the number of degrees of freedom [ default: 1 ]

        prefix (str): The name of the directory containing the q, p, a, and s trajectory data.
    """

#Collect simulation parameters from dyn_params dict...
    nsteps = dyn_params["nsteps"]
    states = dyn_params["states"]
    grid_dims = dyn_params["grid_dims"]
    ndof = dyn_params["ndof"]

    nstates = len(states)
    xmin = plt_params["xmin"]
    xmax = plt_params["xmax"]
    npoints = plt_params["npoints"]

#This is standard Libra convention for wf files...
    data_type1="wfcr"; data_type2="dens"; data_type3="rep_0"

#Determine ntraj from dyn_params grid_dims variable...
    ntraj = 1
    for i in range(len(grid_dims)):
        ntraj *= grid_dims[i]

#Define which timesteps to calculate the wf for...
    lines = plt_params['snaps']

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

    elems = 1
    for i in npoints:
        elems *= i

    b = grid_bounds.flatten()
    wfpts = []
    for i in range(elems):
        index = i
        coords = []

        while index < len(b):
            coords.append(b[index])
            index += elems
        wfpts.append(coords)

#Read the trajectory data and compute the wf on the mesh...
    iline=0
    for line in range(nsteps):
        qdata = qfile.readline().strip().split()
        pdata = pfile.readline().strip().split()
        adata = afile.readline().strip().split()
        coeffs = cfile.readline().strip().split()

        if iline in lines:
            outfile = open(prefix+"/wfc/"+data_type1+"_snap_"+str(line)+"_"+data_type2+"_"+data_type3,"w")

            for pt in wfpts:
                for nn in pt:
                    outfile.write(str(nn)+" ")
                idata = 0
                for state in range(nstates):
                    wf = 0+0j
                    for j in range(ntraj):
                        coeff = complex(float(coeffs[2*(j+state*ntraj)]),
                                        float(coeffs[2*(j+state*ntraj)+1]))

                        gaus = 1.0+0.0j
                        for dof in range(ndof):
                            qj = float(qdata[idata])
                            pj = float(pdata[idata])
                            aj = float(adata[idata])
                            idata +=1

                            val = pt[dof]
                            gaus*=(aj/np.pi)**0.25*np.exp(-aj/2.0* \
                                (val-qj)**2+1j*(pj*(val-qj)))

                        wf+=coeff*gaus
                    outfile.write(str(abs(wf)**2)+" ")
                    #outfile.write(str(np.imag(wf))+" ")
                outfile.write("\n")
        iline+=1

