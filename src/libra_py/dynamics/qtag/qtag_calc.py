"""
..module:: qtag_calc
  :platform: Unix, Windows
  :synopsis: This module contains "ground-level" functions for calculations, mostly output things like energy, norm, etc.

..moduleauthors :: Matthew Dutra
"""

import sys
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
        n (float): Surface population. The imaginary part should be zero.
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


def super_overlap(ndof,ntraj,nstates,qpas):

#    ndof, ntraj = dyn_params['ndof'],['ntraj']
    super_ov = CMATRIX(ntraj,ntraj)

    for i in range(ntraj):
        for j in range(ntraj):
            super_ov.set(i,j,complex(0.0,0.0))

    #Extract the components of the qpas object into their constituent parts: q, p, a, s, surface IDs.
    qvals,pvals,avals,svals=MATRIX(qpas[0]),MATRIX(qpas[1]),MATRIX(qpas[2]),MATRIX(qpas[3])
    surf_ids=qpas[4]

    itot = 0; jtot = 0
    for n1 in range(nstates):

        #Extract the relevent trajectories for the n1-th state...
        traj_on_surf_n1 = [index for index, traj_id in enumerate(surf_ids) if traj_id == n1]
        ntraj_on_surf_n1 = len(traj_on_surf_n1)

        qvals_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        pvals_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        avals_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)
        svals_surf_n1 = MATRIX(ndof,ntraj_on_surf_n1)

        pop_submatrix(qvals,qvals_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(pvals,pvals_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(avals,avals_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)
        pop_submatrix(svals,svals_surf_n1,[dof for dof in range(ndof)],traj_on_surf_n1)

        ii,jj = ntraj_on_surf_n1, ntraj_on_surf_n1
        n12_mat = CMATRIX(ii,jj)

        n12_mat = state_overlap(qvals_surf_n1, pvals_surf_n1, avals_surf_n1, svals_surf_n1, n1, \
                                qvals_surf_n1, pvals_surf_n1, avals_surf_n1, svals_surf_n1, n1)

        for i in range(ii):
            for j in range(jj):
                super_ov.set(i+itot,j+jtot,n12_mat.get(i,j))

        itot += ii
        jtot += jj

    return(super_ov)


def state_overlap(qvals_surf_n1, pvals_surf_n1, avals_surf_n1, svals_surf_n1, n1, \
                  qvals_surf_n2, pvals_surf_n2, avals_surf_n2, svals_surf_n2, n2):

    """Returns the Gaussian overlap matrix *ov_mat*, which stores the complex overlap elements of the two basis functions defined by the rows of the *qpas1* and *qpas2* matrices.

    Args:
        ntraj (integer): Number of trajectories per surface.

        qpas1 (list): List of {q,p,a,s} MATRIX objects for the first set of basis functions.

        qpas2 (list): List of {q,p,a,s} MATRIX objects for the second set of basis functions.

    Returns:
        ov_mat (CMATRIX): The ntraj-by-ntraj complex (Gaussian) overlap matrix of basis functions defined by qpas1 and qpas2.
    """

    ndof = qvals_surf_n1.num_of_rows
    itraj = qvals_surf_n1.num_of_cols
    jtraj = qvals_surf_n2.num_of_cols
    state_ov = CMATRIX(itraj,jtraj)

    for i in range(itraj):
        for j in range(jtraj):
            state_ov.set(i,j,complex(0.0,0.0))

    for i in range(itraj):
#        ii = i*nstates+n1

        qi = qvals_surf_n1.col(i)
        pi = pvals_surf_n1.col(i)
        ai = avals_surf_n1.col(i)
        si = svals_surf_n1.col(i)

        for j in range(i+1):
#            jj = j*nstates+n2

            qj = qvals_surf_n2.col(j)
            pj = pvals_surf_n2.col(j)
            aj = avals_surf_n2.col(j)
            sj = svals_surf_n2.col(j)

            state_ov.set(i,j,gwp_overlap(qi,pi,si,ai/2,qj,pj,sj,aj/2))

            if j != i:
                state_ov.set(j,i,state_ov.get(i,j).conjugate())

    return(state_ov)


def new_old_overlap(ndof,ntraj,nstates,qpaso,qpasn):
    new_old_ov = CMATRIX(ntraj,ntraj)

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
    for n1 in range(nstates):

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

    return(new_old_ov)

def wf_print_1D(ntraj,qpas,c,fileobj):
    """Writes a 1D single-surface wavefunction to a file specified by *fileobj*, computed using the internal function psi(qpas,c,x0) with basis parameters *qpas* and coefficients *c*.

    Args:
        ntraj (integer): Number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

        fileobj (object): An open file object for printing wavefunction data to.

    Returns:
        None.
    """

    x0=MATRIX(1,1)

    for i in np.linspace(-10.0,12.0,num=1000):
        x0.set(0,0,i)
        z=psi(1,ntraj,qpas,c,x0);zstar=np.conj(z)
        print(x0.get(0,0),np.abs(z),(zstar*z).real,z.real,z.imag,sep=' ', end='\n', file=fileobj)

    print(file=fileobj)
    return()

def wf_print_2D(ntraj,qpas,c,fileobj):
    """Writes a 2D single-surface wavefunction to a file specified by *fileobj*, computed using the internal function psi(qpas,c,x0) with basis parameters *qpas* and coefficients *c*.

    Args:
        ntraj (integer): Number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

        fileobj (object): An open file object for printing wavefunction data to.

    Returns:
        None.
    """

    x0=MATRIX(2,1)

    for i in np.linspace(-6.0,6.0,num=100):
        for j in np.linspace(-6.0,6.0,num=100):
            x0.set(0,0,i);x0.set(1,0,j)
            z=psi(2,ntraj,qpas,c,x0);zstar=np.conj(z)
            print(x0.get(0,0),x0.get(1,0),np.abs(z),(zstar*z).real,sep=' ',end='\n', file=fileobj)

    print(file=fileobj)
    return()

def nonad_assemble(chk,m,n,matrix1,matrix2,matrix3):
	"""Returns a matrix *mtot* assembled from *m*-by-*n* component matrices *matrix1*, *matrix2*, and *matrix3*. The dimensions of *mtot* are dependent upon *m* and *n*: if *m*=*n*, then *mtot* is *2m*-by-*2n*; if *m*!=*n*, then *mtot* is *2m*-by-*n*. The variable *chk* designates whether the input (and output) matrices are real or complex.

        Args:
                chk (string): Variable specifying the type of matrix 1, matrix2, and matrix3. Either 'real' or 'cplx'.

                m (integer): Number of rows of input matrices.

                n (integer): Number of columns of input matrices.

                matrix1 (MATRIX/CMATRIX): First of the matrices to be assembled into a total, system matrix.

                matrix2 (MATRIX/CMATRIX): Second of the matrices to be assembled into a total, system matrix.

                matrix3 (MATRIX/CMATRIX): Third of the matrices to be assembled into a total, system matrix.

        Returns:
                mtot (MATRIX/CMATRIX): Output composite matrix.
	"""

	if chk=='real':
		if m==n:
			mtot=MATRIX(2*m,2*n)
		elif m>n:
			mtot=MATRIX(2*m,n)
		elif m<n:
			mtot=MATRIX(m,2*n)

	elif chk=='cplx':
		if m==n:
			mtot=CMATRIX(2*m,2*n)
		elif m>n:
			mtot=CMATRIX(2*m,n)
		elif m<n:
			mtot=CMATRIX(m,2*n)

	if n==m:
		listm=[]; listn=[]
		for i in range(m):
			listm.append(i)
		for j in range(n):
			listn.append(j)

		push_submatrix(mtot,matrix1,listm,listn)
		push_submatrix(mtot,matrix2,[m+lm for lm in listm],[n+ln for ln in listn])
		push_submatrix(mtot,matrix3,listm,[n+ln for ln in listn])
		push_submatrix(mtot,matrix3.T().conj(),[m+lm for lm in listm],listn)

	elif m>n:
		listm=[]; listn=[]
		for i in range(m):
			listm.append(i)
		for j in range(n):
			listn.append(j)

		push_submatrix(mtot,matrix1,listm,listn)
		push_submatrix(mtot,matrix2,[m+lm for lm in listm],listn)

	elif m<n:
		listm=[]; listn=[]
		for i in range(m):
			listm.append(i)
		for j in range(n):
			listn.append(j)

		push_submatrix(mtot,matrix1,listm,listn)
		push_submatrix(mtot,matrix2,listm,[n+ln for ln in listn])

	return(mtot)
