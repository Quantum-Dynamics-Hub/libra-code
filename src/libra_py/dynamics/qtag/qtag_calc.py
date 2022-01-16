"""
..module:: qtag_calc
  :platform: Unix, Windows
  :synopsis: This module contains "ground-level" functions for calculations, mostly output things like energy, norm, etc.

..moduleauthors :: Matthew Dutra
"""

import sys
from liblibra_core import *

import numpy as np

def psi(ndof, ntraj, qpas, c, x0):
	"""Returns the (complex) wavefunction value *wf* at a given point *x0*, calculated using the single-surface basis parameters stored in *qpas* and coefficients *c*.

        Args:
		ndof (integer): Number of degrees of freedom.

		ntraj (integer): Number of trajectories per surface.

                qpas (list): List of {q,p,a,s} MATRIX objects.

                c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

		x0 (MATRIX): The matrix of coordinates [x_1,...,x_ndof] at which the wavefunction value should be calculated.

        Returns:
                wf (complex): Complex value of the wavefunction at x0.
        """

	qvals,pvals,avals,svals=qpas[0],qpas[1],qpas[2],qpas[3]

	wf=0+0j
	for i in range(ntraj):
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

def norm(c,ov):
	"""Returns a single-surface population *n*, calculated from the single-surface basis coefficients *c* and the appropriate overlap matrix *ov* as <c^T|ov|c>. Note that this is equivalent to the norm for a single-surface system.

        Args:
                c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

                ov (CMATRIX): The single-surface overlap matrix for which the population is to be calculated.

        Returns:
                n (float): Surface population. The imaginary part should be zero.
        """

	n=(c.H()*ov*c).get(0).real
	return(n)

def overlap(ntraj,qpas1,qpas2):
	"""Returns the Gaussian overlap matrix *ov_mat*, which stores the complex overlap elements of the two basis functions defined by the rows of the *qpas1* and *qpas2* matrices.

        Args:
		ntraj (integer): Number of trajectories per surface.

                qpas1 (list): List of {q,p,a,s} MATRIX objects for the first set of basis functions.

                qpas2 (list): List of {q,p,a,s} MATRIX objects for the second set of basis functions.

        Returns:
                ov_mat (CMATRIX): The ntraj-by-ntraj complex (Gaussian) overlap matrix of basis functions defined by qpas1 and qpas2.
        """

	ov_mat=CMATRIX(ntraj,ntraj)

	qvals1,pvals1,avals1,svals1=qpas1[0],qpas1[1],qpas1[2],qpas1[3]
	qvals2,pvals2,avals2,svals2=qpas2[0],qpas2[1],qpas2[2],qpas2[3]

	for i in range(ntraj):
		q1,p1,a1,s1=qvals1.col(i),pvals1.col(i),avals1.col(i),svals1.col(i)
		for j in range(ntraj):
			q2,p2,a2,s2=qvals2.col(j),pvals2.col(j),avals2.col(j),svals2.col(j)

			ov_mat.set(i,j,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2/2))
	return(ov_mat)

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
