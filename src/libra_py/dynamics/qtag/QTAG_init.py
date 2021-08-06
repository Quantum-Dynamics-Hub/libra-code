"""
..module:: QTAG_init
  :platform: Unix, Windows
  :synopsis: This module contains functions for initial basis placement.

..moduleauthors :: Matthew Dutra
"""

import sys
import cmath
import math
import os
from liblibra_core import *
import util.libutil as comn

import numpy as np

"""
.. py:function:: grid(ntraj,traj0,wf0)
   Returns the initial basis parameters {q,p,a,s} as an *ntraj*-by-4
   matrix *qpas*, based on the input contained in the dicts *traj0*
   and *wf0*. The placement is evenly spaced across the domain where
   the initial wavefunction has a magnitude greater than *rcut*.
"""

def grid(ntraj,traj0,wf0):
	"""Returns the initial basis parameters {q,p,a,s} as an *ntraj*-by-4 matrix *qpas*, based on the input contained in the dicts *traj0* and *wf0*. The placement is evenly spaced across the domain where the initial wavefunction has a magnitude greater than *rcut*.

        Args:
                ntraj (integer): The number of trajectories per surface.

                traj0 (dictionary): Dictionary containing initialization parameters and keywords.

		wf0 (dictionary): Dictionary containing initial wavepacket conditions.

        Returns:
                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.
        """

	qpas=MATRIX(ntraj,4)
	rcut=traj0['rho']
	a0=traj0['a0']

	xlow=wf0['q']-np.sqrt(-0.5/wf0['a']*np.log(rcut))
	xhi=wf0['q']+np.sqrt(-0.5/wf0['a']*np.log(rcut))
	qs=np.linspace(xlow,xhi,num=ntraj)

	for i in range(ntraj):
		qpas.set(i,0,qs[i])
		qpas.set(i,1,wf0['p'])
		qpas.set(i,2,wf0['a']*a0)
		qpas.set(i,3,0.0)

	return(qpas)

"""
.. py:function:: gaus(ntraj,traj0,wf0)
   Returns the initial basis parameters {q,p,a,s} as an *ntraj*-by-4
   matrix *qpas*, based on the input contained in the dicts *traj0*
   and *wf0*. The placement is randomly chosen from a Gaussian
   distribution centered about the wavepacket maximum with a 
   standard deviation *rho*. The corresponding Gaussian-distributed
   momenta are ordered so that the wavepacket spreads as x increases.
"""

def gaus(ntraj,traj0,wf0):
	"""Returns the initial basis parameters {q,p,a,s} as an *ntraj*-by-4 matrix *qpas*, based on the input contained in the dicts *traj0* and *wf0*. The placement is randomly chosen from a Gaussian distribution centered about the wavepacket maximum with a standard deviation *rho*. The corresponding Gaussian-distributed momenta are ordered so that the wavepacket spreads as x increases.

        Args:
                ntraj (integer): The number of trajectories per surface.

                traj0 (dictionary): Dictionary containing initialization parameters and keywords.

                wf0 (dictionary): Dictionary containing initial wavepacket conditions.

        Returns:
                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

	"""

	qpas=MATRIX(ntraj,4)
	rho=traj0['rho']
	a0=traj0['a0']

	q_gaus=np.sort(np.random.normal(wf0['q'],rho,ntraj))
	p_gaus=np.sort(np.random.normal(wf0['p'],rho,ntraj))

	for i in range(ntraj):
		qpas.set(i,0,q_gaus[i])
		qpas.set(i,1,p_gaus[i])
		qpas.set(i,2,wf0['a']*a0)
		qpas.set(i,3,0.0)

	return(qpas)

"""
.. py:function:: coeffs(ntraj,wf0,qpas,nsurf)
   Returns the projection vector *b* of the initial wavefunction
   with parameters stored in the dict *wf0* onto the basis
   defined by *qpas*. This function assumes the wavefunction is
   located entirely on *nsurf*=1 initially.
"""

def coeffs(ntraj,wf0,qpas,nsurf):
	"""Returns the projection vector *b* of the initial wavefunction with parameters stored in the dict *wf0* onto the basis defined by *qpas*. This function assumes the wavefunction is located entirely on *nsurf*=1 initially.

        Args:
                ntraj (integer): The number of trajectories per surface.

                wf0 (dictionary): Dictionary containing initial wavepacket conditions.

		qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

		nsurf (integer): The number of the desired surface. 1 = ground; 2 = excited.

	Returns:
                b (CMATRIX): Projection vector for the initial wavefunction onto the initial Gaussian basis.
        """

	b=CMATRIX(ntraj,1)

	if nsurf==1:
		for i in range(ntraj):
			qpasi=qpas.row(i)
			q1,p1,a1,s1=qpasi.get(0),qpasi.get(1),qpasi.get(2),qpasi.get(3)
			q2,p2,a2,s2=wf0['q'],wf0['p'],wf0['a'],wf0['s']
		
			b.set(i,0,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2))
	elif nsurf==2:
		for i in range(ntraj):
			b.set(i,0,0+0j)
	return(b)

