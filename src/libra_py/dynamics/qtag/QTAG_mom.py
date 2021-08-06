"""
..module:: QTAG_mom
  :platform: Unix, Windows
  :synopsis: This module contains functions for computing basis momentum.

..moduleauthors :: Matthew Dutra
"""

from liblibra_core import *
from QTAG_config import univ, mom_params
import numpy as np
import QTAG_calc

"""
.. py:function:: momentum(ntraj,qpas,c)
   Returns the momentum *mom* calculated for each basis function 
   according to p=Im(grad(psi)/psi). Also returns the corresponding
   real component *r*, which can be used in updates of the basis
   phase parameter *s*.
"""

def _momentum(ntraj,qpas,c):
	"""Returns the momentum *mom* calculated for each basis function according to p=Im(grad(psi)/psi). Also returns the corresponding real component *r*, which can be used in updates of the basis phase parameter *s*.

        Args:
                ntraj (integer): The number of trajectories per surface.

                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

                c (CMATRIX): The coefficient matrix for the basis.

        Returns:
                mom (MATRIX): The momentum matrix, dimensioned ntraj-by-1.

		r (MATRIX): The complementary real component of the momentum, dimensioned ntraj-by-1.
        """
	
	mom=MATRIX(ntraj,1);r=MATRIX(ntraj,1)

	for i in range(ntraj):
		z=complex(0.0,0.0)
		dz=complex(0.0,0.0)
		qpasi=qpas.row(i)
		q1,p1,a1,s1=qpasi.get(0),qpasi.get(1),qpasi.get(2),qpasi.get(3)
		for j in range(ntraj):
			qpasj=qpas.row(j)
			q2,p2,a2,s2=qpasj.get(0),qpasj.get(1),qpasj.get(2),qpasj.get(3)
			term1=np.exp(-0.5*a2*(q1-q2)**2+1.0j*(p2*(q1-q2)+s2))
			z+=c.get(j)*(a2/np.pi)**0.25*term1
			dz-=(a2*(q1-q2)-1.0j*p2)*c.get(j)*(a2/np.pi)**0.25*term1
		ztot=dz/z
		mom.set(i,0,ztot.imag)
		r.set(i,0,ztot.real)
	return(mom,r)

"""
.. py:function:: lin_fit(ntraj,qpas,c,mom_in,r)
   Returns the momentum *mom_out* and it's gradient *gmom* after linear
   fitting via the internal procedure solve_linsys. The points are
   each weighted by the local wavefunction density. The fitted values
   for *r* and its gradient *gr* are also returned.
"""

def _lin_fit(ntraj,qpas,c,mom_in,r):
	"""Returns the momentum *mom_out* and it's gradient *gmom* after linear fitting via the internal procedure solve_linsys. The points are each weighted by the local wavefunction density. The fitted values for *r* and its gradient *gr* are also returned.

        Args:
                ntraj (integer): The number of trajectories per surface.

                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

                c (CMATRIX): The coefficient matrix for the basis.

		mom_in (MATRIX): The input momentum matrix to be fitted.

		r (MATRIX): The real complement to the input momentum matrix to be fitted.

        Returns:
                mom_out (MATRIX): The linear-fitted momentum matrix, dimensioned ntraj-by-1.

                r (MATRIX): The linear-fitted complementary real component of the momentum, dimensioned ntraj-by-1.

                gmom (MATRIX): The fitted momentum gradient matrix, dimensioned ntraj-by-1.

                gr (MATRIX): The fitted gradient of the complementary real component of the momentum, dimensioned ntraj-by-1.
        """

	mom_out=MATRIX(ntraj,1);gmom=MATRIX(ntraj,1);gr=MATRIX(ntraj,1)
	a=MATRIX(2,2);b=MATRIX(2,2);x=MATRIX(2,2)

	for m in range(2):
		bb1=0+0j;bb2=0+0j
		for i in range(ntraj):
			x0=qpas.get(i,0)
			z=QTAG_calc.psi(qpas,c,x0);zstar=np.conj(z)
			bb1+=mom_in.get(i)*qpas.get(i,0)**(m)*(z*zstar).real
			bb2+=r.get(i)*qpas.get(i,0)**(m)*(z*zstar).real
		b.set(m,0,bb1.real);b.set(m,1,bb2.real)

		for n in range(2):
			aa=0+0j
			for i in range(ntraj):
				x0=qpas.get(i,0)
				z=QTAG_calc.psi(qpas,c,x0);zstar=np.conj(z)
				aa+=qpas.get(i,0)**(m+n)*(z*zstar).real
			a.set(m,n,aa.real)

	solve_linsys(a,b,x,mom_params['beta'],200000)
	for i in range(ntraj):
		aa=x.get(0,0)+x.get(1,0)*qpas.get(i,0)
		bb=x.get(0,1)+x.get(1,1)*qpas.get(i,0)
		mom_out.set(i,0,aa);r.set(i,0,bb)
		gmom.set(i,0,x.get(1,0))
		gr.set(i,0,x.get(1,1))
	return(mom_out,r,gmom,gr)

"""
.. py:function:: mom_raw(qpas,c[,*args])
   Returns the raw (i.e. unfitted and unconvoluted) momentum *mom* and
   corresponding real component *r*. The gradient of each (*gmom* and
   *gr*, respectively) are still calculated via a linear fitting
   procedure defined in the function lin_fit, as they are necessary
   for calculations with adaptable width.
"""

def mom_raw(qpas,c,*args):
	"""Returns the raw (i.e. unfitted and unconvoluted) momentum *mom* and corresponding real component *r*. The gradient of each (*gmom* and *gr*, respectively) are still calculated via a linear fitting procedure defined in the function lin_fit, as they are necessary for calculations with adaptable width.

        Args:
                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

                c (CMATRIX): The coefficient matrix for the basis.

        Returns:
                mom (MATRIX): The unaltered momentum matrix, dimensioned ntraj-by-1.

                r (MATRIX): The unaltered complementary real component of the momentum, dimensioned ntraj-by-1.

                gmom (MATRIX): The fitted momentum gradient matrix, dimensioned ntraj-by-1.

                gr (MATRIX): The fitted gradient of the complementary real component of the momentum, dimensioned ntraj-by-1.
        """

	ntraj=univ['ntraj']
	aaa=MATRIX(ntraj,1);bbb=MATRIX(ntraj,1)
	mom,r=_momentum(ntraj,qpas,c)
	aaa,bbb,gmom,gr=_lin_fit(ntraj,qpas,c,mom,r)
	return(mom,r,gmom,gr)

"""
.. py:function:: mom_fit(qpas,c[,*args])
   Returns the basis momentum *mom* and corresponding real component
   *r* after linear fitting via the function lin_fit. The gradients
   of each (*gmom* and *gr*, respectively) are also returned.
"""

#Compute momentum from linear fitting of values.
def mom_fit(qpas,c,*args):
	"""Returns the basis momentum *mom* and corresponding real component *r* after linear fitting via the function lin_fit. The gradients of each (*gmom* and *gr*, respectively) are also returned.

        Args:
                qpas (MATRIX): The ntraj-by-4 basis parameter matrix.

                c (CMATRIX): The coefficient matrix for the basis.

        Returns:
                mom (MATRIX): The linear-fitted momentum matrix, dimensioned ntraj-by-1.

                r (MATRIX): The linear-fitted complementary real component of the momentum, dimensioned ntraj-by-1.

                gmom (MATRIX): The fitted momentum gradient matrix, dimensioned ntraj-by-1.

                gr (MATRIX): The fitted gradient of the complementary real component of the momentum, dimensioned ntraj-by-1.
	"""

	ntraj=univ['ntraj']
	mom,r=_momentum(ntraj,qpas,c)
	mom,r,gmom,gr=_lin_fit(ntraj,qpas,c,mom,r)
	return(mom,r,gmom,gr)

"""
.. py:function:: mom_avg(qpas1,c1,qpas2,c2)
   Returns the basis momenta *mom* and corresponding real component *r*
   after performing a density-weighted average across both surfaces at each
   point. The averaged momenta are then fitted via the lin_fit function.
   The gradients *gmom* and *gr* are also returned.
"""

#Compute momentum as an average of the two surfaces weighted by the wavefunction on each.
def mom_avg(qpas1,c1,qpas2,c2):
	"""Returns the basis momenta *mom* and corresponding real component *r* after performing a density-weighted average across both surfaces at each point. The averaged momenta are then fitted via the lin_fit function. The gradients *gmom* and *gr* are also returned.

        Args:
                qpas1 (MATRIX): The ntraj-by-4 basis parameter matrix for the functions on surface 1.

                c1 (CMATRIX): The coefficient matrix for the basis on surface 1.

                qpas2 (MATRIX): The ntraj-by-4 basis parameter matrix for the functions on surface 2.

                c2 (CMATRIX): The coefficient matrix for the basis on surface 2.

        Returns:
                mom (MATRIX): The linear-fitted momentum matrix averaged across both surfaces, dimensioned ntraj-by-1.

                r (MATRIX): The linear-fitted complementary real component of the average momentum, dimensioned ntraj-by-1.

                gmom (MATRIX): The fitted momentum gradient matrix, dimensioned ntraj-by-1.

                gr (MATRIX): The fitted gradient of the complementary real component of the momentum, dimensioned ntraj-by-1.
        """

	ntraj=univ['ntraj']
	aaa=MATRIX(ntraj,1);bbb=MATRIX(ntraj,1);mom=MATRIX(ntraj,1);r=MATRIX(ntraj,1)
	mom1,r1=_momentum(ntraj,qpas1,c1)
	mom2,r2=_momentum(ntraj,qpas2,c2)

	for i in range(ntraj):
		z=QTAG_calc.psi(qpas1,c1,qpas1.get(i,0));zs=np.conj(z)
		d1=(z*zs).real
		z=QTAG_calc.psi(qpas2,c2,qpas2.get(i,0));zs=np.conj(z)
		d2=(z*zs).real
		d12=d1+d2
		mom.set(i,0,d1*mom1.get(i,0)/d12+d2*mom2.get(i,0)/d12)
		r.set(i,0,d1*r1.get(i,0)/d12+d2*r2.get(i,0)/d12)

	mom,r,gmom,gr=_lin_fit(ntraj,qpas1,c1,mom,r)
	return(mom,r,gmom,gr)
