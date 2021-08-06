"""
..module:: QTAG_ham
  :platform: Unix, Windows
  :synopsis: This module contains functions useful for computing Hamiltonian matrix elements and the subsequent basis coefficient updates.

..moduleauthors :: Matthew Dutra
"""

import sys
import cmath
import math
import os
from liblibra_core import *
import util.libutil as comn

import numpy as np

from QTAG_config import univ, model_params
import QTAG_calc

"""
.. py:function:: LHA(qpasi,qpasj,nsurf,pot[,**model_params])
   Returns the (complex) value for the potential *v* on an energetic
   surface specified by *nsurf* from two basis functions defined by 
   their parameters *qpasi* and *qpasj*, respectively. The computation
   employs the Local Harmonic Approximation (LHA), which requires a 
   potential function *pot* as well as its first and second derivatives.
   The potential parameters are stored in the dict 'model_params' as
   defined in QTAG_config.
"""

def LHA(qpasi,qpasj,nsurf,pot,**model_params):
	"""Returns the (complex) value for the potential *v* on an energetic surface specified by *nsurf* from two basis functions defined by their parameters *qpasi* and *qpasj*, respectively. The computation employs the Local Harmonic Approximation (LHA), which requires a potential function *pot* as well as its first and second derivatives. The potential parameters are stored in the dict 'model_params' as defined in QTAG_config.

        Args:
                qpasi (MATRIX): The 1-by-4 parameter matrix of the i-th basis function.

                qpasj (MATRIX): The 1-by-4 parameter matrix of the j-th basis function.

		nsurf (integer): Integer specifying the potential surface to be calculated. 1 = ground; 2 = excited; 3 = coupling.

		pot (function object): Function object containing the potential model for the ground and excited states.

        Returns:
                v (complex): The complex potential v computed via the local harmonic approximation between basis functions i and j.
        """

	v=complex(0.0,0.0)
	q1,p1,a1=qpasi.get(0), qpasi.get(1), qpasi.get(2)
	q2,p2,a2=qpasj.get(0), qpasj.get(1), qpasj.get(2)

	z=complex(a1*q1+a2*q2,p2-p1)/(a1+a2)
	xs=[q1,q2]
	for x in xs:
		vx,dvx,d2vx=pot(x,nsurf,**model_params)
		vv0=vx-dvx*x+d2vx/2.0*x**2
		vv1=-d2vx*x+dvx
		vv2=d2vx/2.0
		v+=0.5*(vv0+vv1*z+vv2*(z**2+1.0/(a1+a2)))

	return(v)

"""
.. py:function:: hamiltonian(qpas,ov,nsurf,pot)
   Calculates the single-surface Hamiltonian matrix elements H_ij=<gi|KE+V|gj>,
   computed using the basis parameters stored in *qpas*. This requires the 
   single-surface overlap matrix *ov* as well as the potential function *pot*,
   which is specified in QTAG_config and selected by QTAG_assembler. The surface
   is designated by *nsurf*. Returns the single-surface Hamiltonian *H*.
"""

def hamiltonian(qpas,ov,nsurf,pot):
	"""Calculates the single-surface Hamiltonian matrix elements H_ij=<gi|KE+V|gj>, computed using the basis parameters stored in *qpas*. This requires the single-surface overlap matrix *ov* as well as the potential function *pot*, which is specified in QTAG_config and selected by QTAG_assembler. The surface is designated by *nsurf*. Returns the single-surface Hamiltonian *H*.

        Args:
                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

                ov (CMATRIX): The single-surface ntraj-by-ntraj complex overlap matrix.

                nsurf (integer): Integer specifying the potential surface to be calculated. 1 = ground; 2 = excited; 3 = coupling.

                pot (function object): Function object containing the potential model for the ground and excited states.

        Returns:
                H (CMATRIX): The single-surface ntraj-by-ntraj Hamiltonian matrix.
        """

	ntraj=univ['ntraj']
	H = CMATRIX(ntraj,ntraj)
	q1=MATRIX(1,1); q2=MATRIX(1,1)
	p1=MATRIX(1,1); p2=MATRIX(1,1)
	a1=MATRIX(1,1); a2=MATRIX(1,1)
	s1=MATRIX(1,1); s2=MATRIX(1,1)

	for i in range(ntraj):
		for j in range(ntraj):
			qpasi=qpas.row(i)
			qpasj=qpas.row(j)
			q1.set(0,0,qpasi.get(0));p1.set(0,0,qpasi.get(1));a1.set(0,0,qpasi.get(2));s1.set(0,0,qpasi.get(3))
			q2.set(0,0,qpasj.get(0));p2.set(0,0,qpasj.get(1));a2.set(0,0,qpasj.get(2));s2.set(0,0,qpasj.get(3))
			ke=(-1.0**2/(2.0*univ['mass']))*gwp_kinetic(q1,p1,s1,a1/2,q2,p2,s2,a2/2)

			v=LHA(qpasi,qpasj,nsurf,pot,**model_params)
			H.set(i,j,ke+v*ov.get(i,j))

	return(H)

"""
.. py:function:: coupling(qpas1,qpas2,ov12,nsurf,cplg,pot)
   Calculates the off-diagonal coupling matrix elements V_ij=<gi|V_cpl|gj>,
   where gi and gj are basis functions on separate surfaces defined by *qpas1*
   and *qpas2*, respectively. This requires the overlap matrix *ov12* as well as 
   the potential function *pot* and coupling type *cplg*. The value of *nsurf*
   is set to 3 in the function call to specify the coupling function and its
   derivatives should be used if invoking the LHA (rather than the single
   surfaces 1 or 2). Returns the off-diagonal matrix *H*.
"""

def coupling(qpas1,qpas2,ov12,nsurf,cplg,pot):
	"""Calculates the off-diagonal coupling matrix elements V_ij=<gi|V_cpl|gj>, where gi and gj are basis functions on separate surfaces defined by *qpas1* and *qpas2*, respectively. This requires the overlap matrix *ov12* as well as the potential function *pot* and coupling type *cplg*. The value of *nsurf* is set to 3 in the function call to specify the coupling function and its derivatives should be used if invoking the LHA (rather than the single surfaces 1 or 2). Returns the off-diagonal matrix *H*.

        Args:
                qpas1 (MATRIX): The ntraj-by-4 basis parameter matrix on surface 1.

		qpas2 (MATRIX): The ntraj-by-4 basis parameter matrix on surface 2.

                ov12 (CMATRIX): The dual-surface ntraj-by-ntraj complex overlap matrix.

                nsurf (integer): Integer specifying the potential surface to be calculated. 1 = ground; 2 = excited; 3 = coupling (always set to 3 here).

		cplg (function object): Function object containing the coupling potential.

                pot (function object): Function object containing the potential model for the ground and excited states.

        Returns:
                H (CMATRIX): The coupling-surface ntraj-by-ntraj Hamiltonian matrix.
        """
	ntraj=univ['ntraj']
	H = CMATRIX(ntraj,ntraj)

	for i in range(ntraj):
		for j in range(ntraj):
			qpasi=qpas1.row(i)
			qpasj=qpas2.row(j)

			v=cplg(qpasi,qpasj,nsurf,pot,**model_params)
			H.set(i,j,v*ov12.get(i,j))

	return(H)

"""
.. py:function:: basis_diag(m,H,ov,b)
   Returns the updated basis coefficients for both surfaces, stored in a single
   vector *c_new*, computed as c_new=Z*exp(-i*dt*eps)*Z_dagger*b. The variables
   eps and Z are the eigenvalues and eigenvectors obtained from diagonalizing the
   full *m*-by-*m* Hamiltonian matrix *H* using the solve_eigen internal function.
   Note that the projection vector *b* and the full overlap matrix *ov* are
   required.
"""

def basis_diag(m,H,ov,b):
	"""Returns the updated basis coefficients for both surfaces, stored in a single vector *c_new*, computed as c_new=Z*exp(-i*dt*eps)*Z_dagger*b. The variables eps and Z are the eigenvalues and eigenvectors obtained from diagonalizing the full *m*-by-*m* Hamiltonian matrix *H* using the solve_eigen internal function. Note that the projection vector *b* and the full overlap matrix *ov* are required.

        Args:
                m (integer): Dimension of the eigenvalue and eigenvector matrices. Usually equal to 2*ntraj.

                H (CMATRIX): Total system Hamiltonian, dimension 2*ntraj-by-2*ntraj.

                ov (CMATRIX): Total overlap matrix, dimension 2*ntraj-by-2*ntraj.

                b (CMATRIX): Projection vector, dimensioned 2*ntraj-by-1.

        Returns:
                c_new (CMATRIX): Updated basis coefficient vector for both surfaces. Computed as c_new=Z*exp(-i*dt*eps)*Z_dagger*b.
        """

	dt=univ['dt']

	evals=CMATRIX(m,m)
	evecs=CMATRIX(m,m)
	solve_eigen(H,ov,evals,evecs,0)
	ct=evecs.T().conj()
	c_new=evecs*(exp_(evals,-dt*1.0j))*ct*b

	return(c_new)

"""
.. py:function:: nonad_assemble(chk,m,n,matrix1,matrix2,matrix3)
   Returns a matrix *mtot* assembled from *m*-by-*n* component matrices 
   *matrix1*, *matrix2*, and *matrix3*. The dimensions of *mtot* are
   dependent upon *m* and *n*: if *m*=*n*, then *mtot* is *2m*-by-*2n*; if
   *m*!=*n*, then *mtot* is *2m*-by-*n*. The variable *chk* designates 
   whether the input (and output) matrices are real or complex. 
"""

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
		else:
			mtot=MATRIX(2*m,n)
	elif chk=='cplx':
		if m==n:
			mtot=CMATRIX(2*m,2*n)
		else:
			mtot=CMATRIX(2*m,n)

	listm=[]; listn=[]
	for i in range(m):
		listm.append(i)
	for j in range(n):
		listn.append(j)

	if n==m:
		listm=[]; listn=[]
		for i in range(m):
			listm.append(i)
		for j in range(n):
			listn.append(j)

		push_submatrix(mtot,matrix1,listm,listn)
		push_submatrix(mtot,matrix2,[m+lm for lm in listm],[n+ln for ln in listn])
		push_submatrix(mtot,matrix3,listm,[n+ln for ln in listn])
		push_submatrix(mtot,matrix3,[m+lm for lm in listm],listn)

	else:
		listm=[]; listn=[]
		for i in range(m):
			listm.append(i)
		for j in range(n):
			listn.append(j)

		push_submatrix(mtot,matrix1,listm,listn)
		push_submatrix(mtot,matrix2,[m+lm for lm in listm],listn)

	return(mtot)

