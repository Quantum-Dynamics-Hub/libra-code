"""
..module:: QTAG_calc
  :platform: Unix, Windows
  :synopsis: This module contains "ground-level" functions for calculations, mostly output things like energy, norm, etc.

..moduleauthors :: Matthew Dutra
"""

import sys
import cmath
import math
import os
from liblibra_core import *
import util.libutil as comn

import numpy as np
from QTAG_config import univ

"""
.. py:function:: psi(qpas,c,x0)
   Returns the (complex) wavefunction value *wf* at a given point *x0*,
   calculated using the single-surface basis parameters stored in *qpas* 
   and coefficients *c*.
"""

def psi(qpas,c,x0):
	"""Returns the (complex) wavefunction value *wf* at a given point *x0*, calculated using the single-surface basis parameters stored in *qpas* and coefficients *c*.

        Args:
                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

                c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

		x0 (float): The position at which the wavefunction value should be calculated.

        Returns:
                wf (complex): Complex value of the wavefunction at x0.
        """

	ntraj=univ['ntraj']
	wf=0+0j
	for i in range(ntraj):
		q1,p1,a1,s1=qpas.get(i,0),qpas.get(i,1),qpas.get(i,2), qpas.get(i,3)
		wf+=c.get(i)*(a1/np.pi)**0.25*np.exp(-a1/2.0*(x0-q1)**2+1j*(p1*(x0-q1)+s1))
	return(wf)

"""
.. py:function:: energy(c,H)
   Returns the system energy *e*, calculated from the total basis coefficients
   *c* and the total Hamiltonian *H* as <c^T|H|c>.
"""

def energy(c,H):
	"""Returns the system energy *e*, calculated from the total basis coefficients *c* and the total Hamiltonian *H* as <c^T|H|c>.

        Args:
                c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

                H (CMATRIX): The full system Hamiltonian (both surfaces + coupling).

        Returns:
                e (complex): Energy. The imaginary part should be zero.
        """

	e=c.T().conj()*H*c
	return(e.get(0))

"""
.. py:function:: norm(c,ov)
   Returns a single-surface population *n*, calculated from the single-surface basis coefficients
   *c* and the appropriate overlap matrix *ov* as <c^T|ov|c>. Note that this is equivalent to
   the norm for a single-surface system.
"""

def norm(c,ov):
	"""Returns a single-surface population *n*, calculated from the single-surface basis coefficients *c* and the appropriate overlap matrix *ov* as <c^T|ov|c>. Note that this is equivalent to the norm for a single-surface system.

        Args:
                c (CMATRIX): The ntraj-by-1 complex matrix of basis coefficients.

                ov (CMATRIX): The single-surface overlap matrix for which the population is to be calculated.

        Returns:
                n (complex): Surface population. The imaginary part should be zero.
        """

	n=c.T().conj()*ov*c
	return(n.get(0))

"""
.. py:function:: update(ndim,a,b)
   Returns an updated matrix *a* from an input matrix *b*. The assumed shape of *b* is
   *ntraj*-by-*ndim*.
"""

def update(ndim,a,b):
	"""Returns an updated matrix *a* from an input matrix *b*. The assumed shape of *b* is *ntraj*-by-*ndim*.

        Args:
                ndim (integer): Dimensionality (number of columns) of the matrix being updated. The number of rows is assumed to be ntraj.

                b (MATRIX/CMATRIX): Input MATRIX/CMATRIX.

        Returns:
                a (MATRIX/CMATRIX): Output MATRIX/CMATRIX.
        """

	ntraj=univ['ntraj']
	for n in range(ndim):
		for i in range(ntraj):
			a.set(i,n,b.get(i,n))
	return(a)

"""
.. py:function:: overlap(qpas1,qpas2)
   Returns the Gaussian overlap matrix *ov_mat*, which stores the complex overlap elements of
   the two basis functions defined by the rows of the *qpas1* and *qpas2* matrices.
"""

def overlap(qpas1,qpas2):
	"""Returns the Gaussian overlap matrix *ov_mat*, which stores the complex overlap elements of the two basis functions defined by the rows of the *qpas1* and *qpas2* matrices.

        Args:
                qpas1 (MATRIX): The ntraj-by-4 matrix of basis parameters for the first set of basis functions.

                qpas2 (MATRIX): The ntraj-by-4 matrix of basis parameters for the second set of basis functions.

        Returns:
                ov_mat (CMATRIX): The ntraj-by-ntraj complex (Gaussian) overlap matrix of basis functions defined by qpas1 and qpas2.
        """

	ntraj=univ['ntraj']
	ov_mat=CMATRIX(ntraj,ntraj)

	for i in range(ntraj):
		for j in range(ntraj):
			qpasi=qpas1.row(i)
			qpasj=qpas2.row(j)
			q1,p1,a1,s1=qpasi.get(0),qpasi.get(1),qpasi.get(2), qpasi.get(3)
			q2,p2,a2,s2=qpasj.get(0),qpasj.get(1),qpasj.get(2), qpasj.get(3)

			ov_mat.set(i,j,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2/2))
	return(ov_mat)

"""
.. py:function:: wf_print(qpas,c,fileobj)
   Writes a single-surface wavefunction to a file specified by *fileobj*, computed using the
   internal function psi(qpas,c,x0) with basis parameters *qpas* and coefficients *c*.
"""

def wf_print(qpas,c,fileobj):
	"""Writes a single-surface wavefunction to a file specified by *fileobj*, computed using the internal function psi(qpas,c,x0) with basis parameters *qpas* and coefficients *c*.

        Args:
                qpas (MATRIX): The ntraj-by-1 complex matrix of basis coefficients.

                c (CMATRIX): The full system Hamiltonian (both surfaces + coupling).

		fileobj (object): An open file object for printing wavefunction data to.

	Returns:
		None.
        """

	ntraj=univ['ntraj']

	for x0 in np.linspace(-10.0,12.0,num=1000):
		z=psi(qpas,c,x0);zstar=np.conj(z)
		print(x0,np.abs(z),(zstar*z).real,z.real,z.imag,sep=' ', end='\n', file=fileobj)

	print(file=fileobj)
	return()

