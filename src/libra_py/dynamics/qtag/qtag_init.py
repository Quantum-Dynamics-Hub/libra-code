"""
..module:: qtag_init
  :platform: Unix, Windows
  :synopsis: This module contains functions for initial basis placement.

..moduleauthors :: Matthew Dutra
"""

import sys
import os
from liblibra_core import *

import numpy as np

def qtag_checks(univ,wf0,traj0,mss,model,model_params):
	"""Runs checks to ensure the validity of the user input.

        Args:
                univ (dictionary): Dictionary containing system parameters.

                wf0 (dictionary): Dictionary containing initial wavepacket conditions.

                traj0 (dictionary): Dictionary containing initialization parameters and keywords.

		mss (dictionary): Dictionary containing multi-surface scheme keywords.

                model (dictionary): Dictionary containing containing keywords for the potential energy model.

		model_params (dictionary): Dictionary containing the potential energy parameters.
        Returns:
		None.
	"""

	ntraj,ndof=univ['ntraj'],univ['ndof']
	if len(univ['mass']) != ndof:
		sys.exit("List of masses does not match given number of DoFs!")
	if len(wf0['q']) != ndof:
		sys.exit("Initial WF q's do not match given number of DoFs!")
	if len(wf0['p']) != ndof:
		sys.exit("Initial WF p's do not match given number of DoFs!")
	if len(wf0['a']) != ndof:
		sys.exit("Initial WF a's do not match given number of DoFs!")
	if len(wf0['s']) != ndof:
		sys.exit("Initial WF s's do not match given number of DoFs!")
	if len(traj0['grid_dims']) != ndof and traj0['placement']=='grid':
		sys.exit("List of grid_dims does not match given number of DoFs!")
	if len(traj0['a0']) != ndof:
		sys.exit("List of a0 values does not match given number of DoFs!")
	if model['rep'] == 'diabatic':
		model['rep'] = 'diab'
	if model['rep'] == 'adiabatic':
		model['rep'] = 'adiab'
	if model['coupling'] == 'exact':
		model['coupling'] = 'exact_gauss_cpl'
	if model['coupling'] == 'none':
		model['coupling'] = 'no_coupling'
	if model['pot_type'] == 'HO':
		for i in range(ndof):
			model_params['d2'][i]=model_params['d2'][i]*np.sqrt(model_params['k1'][i])
	if np.prod(traj0['grid_dims'])!=ntraj:
		sys.exit("Error in grid_dims: number does not match given value for ntraj!")

	return()	

def grid(ndof,ntraj,traj0,wf0):
	"""Returns the initial basis parameters {q,p,a,s} as an *ntraj*-by-4 matrix *qpas*, based on the input contained in the dicts *traj0* and *wf0*. The placement is evenly spaced across the domain where the initial wavefunction has a magnitude greater than *rcut*.

        Args:
		ndof (integer): The number of degrees of freedom.

                ntraj (integer): The number of trajectories per surface.

                traj0 (dictionary): Dictionary containing initialization parameters and keywords.

		wf0 (dictionary): Dictionary containing initial wavepacket conditions.

        Returns:
                qpas (list): List of {q,p,a,s} MATRIX objects.
        """

	qvals=MATRIX(ndof,ntraj)
	pvals=MATRIX(ndof,ntraj)
	avals=MATRIX(ndof,ntraj)
	svals=MATRIX(ndof,ntraj)

	rcut=traj0['rho']
	a0=traj0['a0']
	grid_dims=traj0['grid_dims']

	qlo,qhi=[],[]
	for i in range(ndof):
        	xlow=wf0['q'][i]-np.sqrt(-0.5/wf0['a'][i]*np.log(rcut))
        	xhi=wf0['q'][i]+np.sqrt(-0.5/wf0['a'][i]*np.log(rcut))
        	qlo.append(xlow)
        	qhi.append(xhi)

	grid=np.mgrid[tuple(slice(qlo[i],qhi[i],complex(0,grid_dims[i])) for i in range(ndof))]

	elems=1
	for i in grid_dims:
        	elems*=i

	b=grid.flatten()
	qs=[]
	for i in range(elems):
        	index=i
        	coords=[]
        	while index < len(b):
                	coords.append(b[index])
                	index+=elems
        	qs.append(coords)

	for i in range(ndof):
		for j in range(ntraj):
			qvals.set(i,j,qs[j][i])
			pvals.set(i,j,wf0['p'][i])
			avals.set(i,j,wf0['a'][i]*a0[i])
			svals.set(i,j,0.0)

	qpas=[qvals,pvals,avals,svals]
	return(qpas)

def gaussian(ndof,ntraj,traj0,wf0):
	"""Returns the initial basis parameters {q,p,a,s} as an *ntraj*-by-4 matrix *qpas*, based on the input contained in the dicts *traj0* and *wf0*. The placement is randomly chosen from a Gaussian distribution centered about the wavepacket maximum with a standard deviation *rho*. The corresponding Gaussian-distributed momenta are ordered so that the wavepacket spreads as x increases.

        Args:
		ndof (integer): The number of degrees of freedom.

                ntraj (integer): The number of trajectories per surface.

                traj0 (dictionary): Dictionary containing initialization parameters and keywords.

                wf0 (dictionary): Dictionary containing initial wavepacket conditions.

        Returns:
                qpas (list): List of {q,p,a,s} MATRIX objects.

	"""

	qvals=MATRIX(ndof,ntraj)
	pvals=MATRIX(ndof,ntraj)
	avals=MATRIX(ndof,ntraj)
	svals=MATRIX(ndof,ntraj)

	rho=traj0['rho']
	a0=traj0['a0']

	q_gaus=np.sort(np.random.normal(wf0['q'],rho,ntraj))
	p_gaus=np.sort(np.random.normal(wf0['p'],rho,ntraj))

	for i in range(ndof):
		for j in range(ntraj):
			qvals.set(i,j,q_gaus[i])
			pvals.set(i,j,p_gaus[i])
			avals.set(i,j,wf0['a']*a0)
			svals.set(i,j,0.0)

	qpas=[qvals,pvals,avals,svals]
	return(qpas)

def coeffs(ndof,ntraj,wf0,qpas,nsurf):
	"""Returns the projection vector *b* of the initial wavefunction with parameters stored in the dict *wf0* onto the basis defined by *qpas*. This function assumes the wavefunction is located entirely on *nsurf*=1 initially.

        Args:
		ndof (integer): The number of degrees of freedom.

                ntraj (integer): The number of trajectories per surface.

                wf0 (dictionary): Dictionary containing initial wavepacket conditions.

		qpas (list): List of {q,p,a,s} MATRIX objects.

		nsurf (integer): The number of the desired surface. 1 = ground; 2 = excited.

	Returns:
                b (CMATRIX): Projection vector for the initial wavefunction onto the initial Gaussian basis.
        """

	q2,p2,a2,s2=MATRIX(ndof,1),MATRIX(ndof,1),MATRIX(ndof,1),MATRIX(ndof,1)
	b=CMATRIX(ntraj,1)
	qvals,pvals,avals,svals=qpas[0],qpas[1],qpas[2],qpas[3]

	for i in range(ndof):
		q2.set(i,0,wf0['q'][i])
		p2.set(i,0,wf0['p'][i])
		a2.set(i,0,wf0['a'][i])
		s2.set(i,0,wf0['s'][i])

	if nsurf==1:
		for i in range(ntraj):
			q1,p1,a1,s1=qvals.col(i),pvals.col(i),avals.col(i),svals.col(i)
			b.set(i,0,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2))
	elif nsurf==2:
		for i in range(ntraj):
			b.set(i,0,0+0j)
#			q1,p1,a1,s1=qvals.col(i),pvals.col(i),avals.col(i),svals.col(i)
#			b.set(i,0,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2)/np.sqrt(2.0))
	elif nsurf==3:
		for i in range(ntraj):
			if i%2==0:
				q1,p1,a1,s1=qvals.col(i),pvals.col(i),avals.col(i),svals.col(i)
				b.set(i,0,gwp_overlap(q1,p1,s1,a1/2,q2,p2,s2,a2))
			else:
				b.set(i,0,0+0j)

	return(b)

