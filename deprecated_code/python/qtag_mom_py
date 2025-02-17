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
..module:: qtag_mom
  :platform: Unix, Windows
  :synopsis: This module contains functions for computing basis momentum.

..moduleauthors :: Matthew Dutra, Alexey Akimov
"""

from liblibra_core import *
from libra_py import data_outs
import numpy as np

from . import qtag_calc

import util.libutil as comn

"""
def _momentum(ndof, ntraj_on_surf, qpas, c):
    Returns the momentum *mom* calculated for each basis function according to p=Im(grad(psi)/psi). 
       Also returns the corresponding real component *r*, which can be used in updates of the basis phase parameter *s*.

    Args:

        ndof (integer): The number of degrees of freedom.

        ntraj_on_surf (integer): The number of trajectories per surface.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        c (CMATRIX(ntraj_on_surf x 1) ): The coefficient matrix for the basis.

    Returns:

        mom (MATRIX): The momentum matrix, dimensioned ndof-by-ntraj_on_surf.

        r (MATRIX): The complementary real component of the momentum, dimensioned ndof-by-ntraj.


    mom=MATRIX(ndof,ntraj_on_surf)
    r=MATRIX(ndof,ntraj_on_surf) 
    dz=CMATRIX(ndof,1)

    qvals,pvals,avals,svals=MATRIX(qpas[0]),MATRIX(qpas[1]),MATRIX(qpas[2]),MATRIX(qpas[3])

    for i in range(ntraj_on_surf):
        z=complex(0.0,0.0)
        for n in range(ndof):
            dz.set(n,0,0+0j)

        q1,p1,a1,s1=qvals.col(i),pvals.col(i),avals.col(i),svals.col(i)

        for j in range(ntraj_on_surf):
            term0,term1 = 1.0,1.0
            q2,p2,a2,s2 = qvals.col(j),pvals.col(j),avals.col(j),svals.col(j)
            for n in range(ndof):
                term0 *= (a2.get(n)/np.pi)**0.25
                term1 *= np.exp(-0.5*a2.get(n)*(q1.get(n)-q2.get(n))**2+1.0j*(p2.get(n)*(q1.get(n)-q2.get(n))+s2.get(n)))
            z += c.get(j)*term0*term1

            for n in range(ndof):
                dzt = dz.get(n)-(a2.get(n)*(q1.get(n)-q2.get(n))-1.0j*p2.get(n))*c.get(j)*term0*term1
                dz.set(n,dzt)	

        for n in range(ndof):
            mom.set(n,i,(dz.get(n)/z).imag)
            r.set(n,i,(dz.get(n)/z).real)

    return (mom,r)


def _lin_fitting(ndof,ntraj_on_surf,qvals,qpas,c,mom_in,r,d_weight,beta):
    Returns the momentum *mom_out* and it's gradient *gmom* after linear fitting via the 
       internal procedure solve_linsys. The points are each weighted by the local wavefunction density. 
       The fitted values for *r* and its gradient *gr* are also returned.

    Args:
        ndof (integer): The number of degrees of freedom.

        ntraj_on_surf (integer): The number of trajectories per surface.

		qvals (MATRIX): ndof-by-ntraj_on_surf MATRIX object containing trajectory positions.

                qpas (list): List of {q,p,a,s} MATRIX objects.

                c (CMATRIX(ntraj_on_surf, 1)): The coefficient matrix for the basis.

		mom_in (MATRIX): The input momentum matrix to be fitted.

		r (MATRIX): The real complement to the input momentum matrix to be fitted.

		d_weight (integer): Integer specifying if the fitted momentum should be density-weighted.

		beta (real): Parameter from the univ dictionary specifying tolerance for the linear fitting algorithm.

        Returns:
                mom_out (MATRIX): The linear-fitted momentum matrix, dimensioned ndof-by-ntraj_on_surf.

                r (MATRIX): The linear-fitted complementary real component of the momentum, dimensioned ndof-by-ntraj_on_surf.

                gmom (MATRIX): The fitted momentum gradient matrix, dimensioned ndof-by-ntraj_on_surf.

                gr (MATRIX): The fitted gradient of the complementary real component of the momentum, dimensioned ndof-by-ntraj_on_surf.

    aaa=MATRIX(1,ntraj_on_surf)
    bbb=MATRIX(1,ntraj_on_surf)
    ccc=MATRIX(1,ntraj_on_surf)
    ddd=MATRIX(1,ntraj_on_surf)

    a=MATRIX(2,2);b=MATRIX(2,2);x=MATRIX(2,2)

    for m in range(2):

        bb1=0+0j;bb2=0+0j

        for i in range(ntraj_on_surf):
            x0 = qpas[0].col(i)

            if d_weight == 1:

                #CMATRIX qtag_psi(MATRIX& q, MATRIX& q1, MATRIX& p1, MATRIX& alp1, MATRIX& s1, CMATRIX& Coeff);
                z = qtag_calc.psi(ndof,ntraj_on_surf,qpas,c,x0);
                zstar = np.conj(z)
            else:
                z=1+0j;zstar=1-0j

                bb1+=mom_in.get(i)*qvals.get(i)**(m)*(z*zstar).real
                bb2+=r.get(i)*qvals.get(i)**(m)*(z*zstar).real
        b.set(m,0,bb1.real);b.set(m,1,bb2.real)

        for n in range(2):
            aa=0+0j
            for i in range(ntraj_on_surf):
                x0=qpas[0].col(i)

                if d_weight == 1:
                    z=qtag_calc.psi(ndof,ntraj_on_surf,qpas,c,x0);zstar=np.conj(z)
                else:
                    z=1+0j;zstar=1-0j

                aa+=qvals.get(i)**(m+n)*(z*zstar).real
            a.set(m,n,aa.real)

    solve_linsys(a,b,x,beta,200000)
    for i in range(ntraj_on_surf):
        aa=x.get(0,0)+x.get(1,0)*qvals.get(i)
        bb=x.get(0,1)+x.get(1,1)*qvals.get(i)
        aaa.set(i,aa);bbb.set(i,bb)
        ccc.set(i,x.get(1,0))
        ddd.set(i,x.get(1,1))

    return(aaa,bbb,ccc,ddd) # Im, Re, dIm, dRe

def unmodified(ndof,ntraj_on_surf,beta,qpas,c,*args):
    Returns the raw (i.e. unfitted and unconvoluted) momentum *mom* and corresponding real component *r*. The gradient of each (*gmom* and *gr*, respectively) are still calculated via a linear fitting procedure defined in the function lin_fit, as they are necessary for calculations with adaptable width.

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        beta (real): Parameter from the univ dictionary specifying tolerance for the linear fitting algorithm.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        c (CMATRIX): The coefficient matrix for the basis.

    Returns:
        mom (MATRIX): The unaltered momentum matrix, dimensioned ndof-by-ntraj.

        r (MATRIX): The unaltered complementary real component of the momentum, dimensioned ndof-by-ntraj.

        gmom (MATRIX): The fitted momentum gradient matrix, dimensioned ndof-by-ntraj.

        gr (MATRIX): The fitted gradient of the complementary real component of the momentum, dimensioned ndof-by-ntraj.

    aaa=MATRIX(1,ntraj);bbb=MATRIX(1,ntraj)
    gmom=MATRIX(ndof,ntraj);gr=MATRIX(ndof,ntraj)
    mom,r=_momentum(ndof,ntraj,qpas,c)
    for i in range(ndof):
        qvals=qpas[0].row(i)
        aaa,bbb,ccc,ddd=_lin_fitting(ndof,ntraj,qvals,qpas,c,mom.row(i),r.row(i),1,beta)
        for j in range(ntraj):
            gmom.set(i,j,ccc.get(j));gr.set(i,j,ddd.get(j))

    return(mom,r,gmom,gr)

def lin_fit(ndof,ntraj_on_surf,beta,qpas,c,*args):
    Returns the basis momentum *mom* and corresponding real component *r* after linear fitting via the function lin_fit. The gradients of each (*gmom* and *gr*, respectively) are also returned.

    Args:
        univ (dictionary): Dictionary containing various system parameters.

        beta (real): Parameter from the univ dictionary specifying tolerance for the linear fitting algorithm.

        qpas (list): List of {q,p,a,s} MATRIX objects.

        c (CMATRIX): The coefficient matrix for the basis.

    Returns:
        mom (MATRIX): The linear-fitted momentum matrix, dimensioned ndof-by-ntraj.

        r (MATRIX): The linear-fitted complementary real component of the momentum, dimensioned ndof-by-ntraj.

        gmom (MATRIX): The fitted momentum gradient matrix, dimensioned ndof-by-ntraj.

        gr (MATRIX): The fitted gradient of the complementary real component of the momentum, dimensioned ndof-by-ntraj.
    

    gmom=MATRIX(ndof,ntraj_on_surf);gr=MATRIX(ndof,ntraj_on_surf)
    mom,r=_momentum(ndof,ntraj_on_surf,qpas,c)

    for i in range(ndof):
        qvals=qpas[0].row(i)
        aaa,bbb,ccc,ddd=_lin_fitting(ndof,ntraj_on_surf,qvals,qpas,c,mom.row(i),r.row(i),1,beta)
        for j in range(ntraj_on_surf):
            mom.set(i,j,aaa.get(j));r.set(i,j,bbb.get(j))
            gmom.set(i,j,ccc.get(j));gr.set(i,j,ddd.get(j))

    return(mom,r,gmom,gr)
"""
