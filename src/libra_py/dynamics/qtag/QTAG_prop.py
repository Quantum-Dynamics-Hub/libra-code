"""
..module:: QTAG_pots
  :platform: Unix, Windows
  :synopsis: This module contains functions for basis propagation in different multi-surface schemes (mss).

..moduleauthors :: Matthew Dutra
"""

import sys
import cmath
import math
import os
from liblibra_core import *
from libra_py import data_outs

import numpy as np
import QTAG_calc
from QTAG_config import univ, mss

def sync(mom_calc,props,qpas1,c1_new,qpas2,c2_new,norm2):
	"""Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection vectors *b1* and *b2*, where the motion of both sets of functions are synced to the lower energetic surface while the density on the upper surface is less than a threshold value specified by the *decpl* parameter. Also necessary are the functions for calculating momentum (*mom_calc*) and basis updates (*props*) as selected by the QTAG_assembler module.

        Args:
                mom_calc (function object): The function object needed to calculate the momentum, as defined by QTAG_config.

                props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

		qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

		c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

                qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

                c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

		norm2 (float): The population on surface 2.

        Returns:
                qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

                qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

                b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

                b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
	 """

	ndof,ntraj,dt,mass=univ['ndof'],univ['ntraj'],univ['dt'],univ['mass']
	decpl=mss['decpl']
	qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]
	qvals,pvals,avals,svals=qpas1[0],qpas1[1],qpas1[2],qpas1[3]
	qvalsn,pvalsn,avalsn,svalsn=MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj),MATRIX(ndof,ntraj)

	mom,r,gmom,gr=mom_calc(qpas1,c1_new)

	for i in range(ndof):
		qn=qprop(i,qvals.row(i),mom.row(i));pn=pprop(mom.row(i));an=aprop(i,avals.row(i),gmom.row(i));sn=sprop(i,svals.row(i))
		for j in range(ntraj):
			qvalsn.set(i,j,qn.get(j))
			pvalsn.set(i,j,pn.get(j))
			avalsn.set(i,j,an.get(j))
			svalsn.set(i,j,sn.get(j))

	qpas1n=[qvalsn,pvalsn,avalsn,svalsn]
	ov_no=QTAG_calc.overlap(qpas1n,qpas1)
	b1=ov_no*c1_new

	if norm2 < decpl:
		qpas2n=qpas1n
	else:
		qvals,pvals,avals,svals=qpas2[0],qpas2[1],qpas2[2],qpas2[3]
		mom,r,gmom,gr=mom_calc(qpas2,c2_new)
		for i in range(ndof):
			qn=qprop(i,qvals.row(i),mom.row(i));pn=pprop(mom.row(i));an=aprop(i,avals.row(i),gmom.row(i));sn=sprop(i,svals.row(i))
			for j in range(ntraj):
				qvals.set(i,j,qn.get(i,j))
				pvals.set(i,j,pn.get(i,j))
				avals.set(i,j,an.get(i,j))
				svals.set(i,j,sn.get(i,j))

		qpas2n=[qvals,pvals,avals,svals]

	ov_no=QTAG_calc.overlap(qpas2n,qpas2)
	b2=ov_no*c2_new

	return(qpas1n,qpas2n,b1,b2)

def fixed():
	return()

def mean_field(mom_calc,props,qpas1,c1_new,qpas2,c2_new,*args):
	"""Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection vectors *b1* and *b2*, where the basis parameters are calculated according to an "average" surface, defined by the average momentum calculed in the *mom_avg* function. Also necessary are the functions for calculating momentum (*mom_calc*, although specifically mom_avg here) and basis updates (*props*) as selected by the QTAG_assembler module.

        Args:
                mom_calc (function object): The function object needed to calculate the momentum, as defined by QTAG_config.

                props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

                qpas1 (list): List of {q,p,a,s} MATRIX objects for surface 1.

                c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

                qpas2 (list): List of {q,p,a,s} MATRIX objects for surface 2.

                c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

                norm2 (float): The population on surface 2.

        Returns:
                qpas1n (list): List of updated {q,p,a,s} MATRIX objects for surface 1.

                qpas2n (list): List of updated {q,p,a,s} MATRIX objects for surface 2.

                b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

                b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
        """

	ndof,ntraj,dt,mass=univ['ndof'],univ['ntraj'],univ['dt'],univ['mass']
	decpl=mss['decpl']
	qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]
	qvals,pvals,avals,svals=qpas1[0],qpas1[1],qpas1[2],qpas1[3]

	mom,r,gmom,gr=mom_calc(qpas1,c1_new,qpas2,c2_new)
	qn=qprop(qvals,mom);pn=pprop(mom);an=aprop(avals,gmom);sn=sprop(svals)

	for i in range(ndof):
		for j in range(ntraj):
			qvals.set(i,j,qn.get(i,j))
			pvals.set(i,j,pn.get(i,j))
			avals.set(i,j,an.get(i,j))
			svals.set(i,j,sn.get(i,j))

	qpas1n=[qvals,pvals,avals,svals]; qpas2n=[qvals,pvals,avals,svals]

	ov_no=QTAG_calc.overlap(qpas1n,qpas1)
	b1=ov_no*c1_new

	ov_no=QTAG_calc.overlap(qpas2n,qpas2)
	b2=ov_no*c2_new

	return(qpas1n,qpas2n,b1,b2)