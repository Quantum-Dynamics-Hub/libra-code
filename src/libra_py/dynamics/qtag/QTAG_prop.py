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
import util.libutil as comn

import numpy as np
import QTAG_calc
from QTAG_config import univ, mss

"""
.. py:function:: mss_sync(mom_calc,props,qpas1,c1_new,qpas2,c2_new,norm2)
   Returns the values for the new basis parameter matrices on surfaces 1
   (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection
   vectors *b1* and *b2*, where the motion of both sets of functions are
   synced to the lower energetic surface while the density on the upper
   surface is less than a threshold value specified by the *decpl* parameter.
   Also necessary are the functions for calculating momentum (*mom_calc*)
   and basis updates (*props*) as selected by the QTAG_assembler module.
"""

def mss_sync(mom_calc,props,qpas1,c1_new,qpas2,c2_new,norm2):
	"""Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection vectors *b1* and *b2*, where the motion of both sets of functions are synced to the lower energetic surface while the density on the upper surface is less than a threshold value specified by the *decpl* parameter. Also necessary are the functions for calculating momentum (*mom_calc*) and basis updates (*props*) as selected by the QTAG_assembler module.

        Args:
                mom_calc (function object): The function object needed to calculate the momentum, as defined by QTAG_config.

                props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

		qpas1 (MATRIX): The ntraj-by-4 basis parameter matrix for surface 1.

		c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

                qpas2 (MATRIX): The ntraj-by-4 basis parameter matrix for surface 2.

                c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

		norm2 (float): The population on surface 2.

        Returns:
                qpas1n (MATRIX): The updated ntraj-by-4 basis parameter matrix for surface 1.

                qpas2n (MATRIX): The updated ntraj-by-4 basis parameter matrix for surface 2.

                b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

                b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
	 """

	ntraj,dt,mass=univ['ntraj'],univ['dt'],univ['mass']
	decpl=mss['decpl']
	qpas1n=MATRIX(ntraj,4);qpas2n=MATRIX(ntraj,4)
	qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]

	mom,r,gmom,gr=mom_calc(qpas1,c1_new)
	qn=qprop(qpas1.col(0),mom);pn=pprop(mom);an=aprop(qpas1.col(2),gmom);sn=sprop(qpas1.col(3))
	for i in range(ntraj):
		qpas1n.set(i,0,qn.get(i))
		qpas1n.set(i,1,pn.get(i))
		qpas1n.set(i,2,an.get(i))
		qpas1n.set(i,3,sn.get(i))

	ov_no=QTAG_calc.overlap(qpas1n,qpas1)
	b1=ov_no*c1_new

	if norm2 < decpl:
		qpas2n=qpas1n
	else:
		mom,r,gmom,gr=mom_calc(qpas2,c2_new)
		qn=qprop(qpas2.col(0),mom);pn=pprop(mom);an=aprop(qpas2.col(2),gmom);sn=sprop(qpas2.col(3))
		for i in range(ntraj):
			qpas2n.set(i,0,qn.get(i))
			qpas2n.set(i,1,pn.get(i))
			qpas2n.set(i,2,an.get(i))
			qpas2n.set(i,3,sn.get(i))

	ov_no=QTAG_calc.overlap(qpas2n,qpas2)
	b2=ov_no*c2_new

	return(qpas1n,qpas2n,b1,b2)

"""
.. py:function:: mss_fixed(mom_calc,props,qpas1,c1_new,qpas2,c2_new,norm2)
   Returns the values for the new basis parameter matrices on surfaces 1
   (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection
   vectors *b1* and *b2*, where the basis on the upper surface is fixed
   in the coupling region until there is sufficient wavefunction density
   on the excited surface, as defined by the parameter *decpl*. Also 
   necessary are the functions for calculating momentum (*mom_calc*)
   and basis updates (*props*) as selected by the QTAG_assembler module.
   Note that this method is not yet implemented.
"""

def mss_fixed():
	return()

"""
.. py:function:: mss_sync(mom_calc,props,qpas1,c1_new,qpas2,c2_new[,*args])
   Returns the values for the new basis parameter matrices on surfaces 1
   (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection
   vectors *b1* and *b2*, where the basis parameters are calculated according
   to an "average" surface, defined by the average momentum calculed in the
   *mom_avg* function. Also necessary are the functions for calculating 
   momentum (*mom_calc*, although specifically mom_avg here) and basis updates 
   (*props*) as selected by the QTAG_assembler module.
"""

def mss_mf(mom_calc,props,qpas1,c1_new,qpas2,c2_new,*args):
	"""Returns the values for the new basis parameter matrices on surfaces 1 (*qpas1n*) and 2 (*qpas2n*), as well as their corresponding projection vectors *b1* and *b2*, where the basis parameters are calculated according to an "average" surface, defined by the average momentum calculed in the *mom_avg* function. Also necessary are the functions for calculating momentum (*mom_calc*, although specifically mom_avg here) and basis updates (*props*) as selected by the QTAG_assembler module.

        Args:
                mom_calc (function object): The function object needed to calculate the momentum, as defined by QTAG_config.

                props (list): A list of function objects for propagating basis parameters {q,p,a,s}.

                qpas1 (MATRIX): The ntraj-by-4 basis parameter matrix for surface 1.

                c1_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 1.

                qpas2 (MATRIX): The ntraj-by-4 basis parameter matrix for surface 2.

                c2_new (CMATRIX): The ntraj-by-1 complex coefficient matrix for the basis on surface 2.

                norm2 (float): The population on surface 2.

        Returns:
                qpas1n (MATRIX): The updated ntraj-by-4 basis parameter matrix for surface 1.

                qpas2n (MATRIX): The updated ntraj-by-4 basis parameter matrix for surface 2.

                b1 (CMATRIX): The updated projection vector of the old basis (defined by qpas1) onto the new basis (defined by qpas1n).

                b2 (CMATRIX): The updated projection vector of the old basis (defined by qpas2) onto the new basis (defined by qpas2n).
        """

	ntraj,dt,mass=univ['ntraj'],univ['dt'],univ['mass']
	decpl=mss['decpl']
	qpas1n=MATRIX(ntraj,4);qpas2n=MATRIX(ntraj,4)
	qprop=props[0];pprop=props[1];aprop=props[2];sprop=props[3]

	mom,r,gmom,gr=mom_calc(qpas1,c1_new,qpas2,c2_new)
	qn=qprop(qpas1.col(0),mom);pn=pprop(mom);an=aprop(qpas1.col(2),gmom);sn=sprop(qpas1.col(3))
	for i in range(ntraj):
		qpas1n.set(i,0,qn.get(i)); qpas2n.set(i,0,qn.get(i))
		qpas1n.set(i,1,pn.get(i)); qpas2n.set(i,1,pn.get(i))
		qpas1n.set(i,2,an.get(i)); qpas2n.set(i,2,an.get(i))
		qpas1n.set(i,3,sn.get(i)); qpas2n.set(i,3,sn.get(i))

	ov_no=QTAG_calc.overlap(qpas1n,qpas1)
	b1=ov_no*c1_new

	ov_no=QTAG_calc.overlap(qpas2n,qpas2)
	b2=ov_no*c2_new

	return(qpas1n,qpas2n,b1,b2)
