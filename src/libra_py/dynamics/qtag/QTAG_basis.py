"""
..module:: QTAG_basis
  :platform: Unix, Windows
  :synopsis: This module contains update functions for the {q,p,a,s} basis components.

..moduleauthors :: Matthew Dutra
"""

from liblibra_core import *
from QTAG_config import univ

def frozen(i,param_in,*args):
	"""The non-update function.
        Args:
		i (integer) : The DoF being considered

                param_in (MATRIX): 1-by-ntraj basis parameter matrix for q, p, a, or s.

	Returns:
		param_in (MATRIX): The same input parameter matrix q, p, a, or s.
        """

	return(param_in)

def q_update(ndim,q_in,mom_in):
	"""Returns the updated position *q_out* from an input position *q_in* and the momentum *mom*, calculated via q_out = q_in+mom*dt/m (i.e. symplectic).

        Args:
		ndim (integer): The dimension along which q is being updated.

                q_in (MATRIX): The 1-by-ntraj basis position MATRIX q to be updated.

		mom_in (MATRIX): The 1-by-ntraj momentum MATRIX p corresponding to q.
        Returns:
                q_out (MATRIX): Updated 1-by-ntraj basis position MATRIX q.
        """

	ntraj=univ['ntraj']
	q_out=MATRIX(1,ntraj)

	for i in range(ntraj):
		qdt = q_in.get(i)+mom_in.get(i)*univ['dt']/univ['mass'][ndim]
		q_out.set(i,qdt)
	return(q_out)

def p_update(mom):
	"""Returns the updated basis parameter *p_out* from the momentum *mom*. Currently these two things are equal, but they don't necessarily need to be.

        Args:
                mom (MATRIX): The 1-by-ntraj basis momentum MATRIX.

        Returns:
                p_out (MATRIX): Updated 1-by-ntraj basis parameter MATRIX p.
        """

	p_out = 1.0*mom
	return(p_out)

def a_update(ndim,a_in,gmom):
	"""Returns the updated basis width *a_out*, calculated from the old basis width *a_in* and the momentum gradient *gmom*.

        Args:
		ndim (integer): The dimension along which a is being updated.

                a_in (MATRIX): The 1-by-ntraj basis width MATRIX a to be updated.

		gmom (MATRIX): The 1-by-ntraj momentum gradient MATRIX storing the values of (dp/dx) at the basis positions corresponding to MATRIX a.
        Returns:
                a_out (MATRIX): Updated 1-by-ntraj basis width MATRIX a.
        """

	ntraj=univ['ntraj']
	a_out,an1 = MATRIX(1,ntraj),MATRIX(1,ntraj)
	an1.dot_product(a_in,gmom)
	for i in range(ntraj):
		adt=a_in.get(i)-2.0*an1.get(i)*univ['dt']/univ['mass'][ndim]
		a_out.set(i,adt)
	return(a_out)

def param_check(basis):
	"""Checks the type of propagation requested for each basis parameter {q,p,a,s} from the basis dictionary.

        Args:
                basis (dictionary): The input dictionary containing keywords for the basis parameter propagation.
        Returns:
                props (list of function objects): List containing the functions for {q,p,a,s} parameter updates.
	"""

	if basis['qtype'] == 'adpt':
		qprop=q_update
	elif basis['qtype'] == 'frzn':
		qprop=frozen
	else:
		sys.exit("Unrecognized keyword in basis qtype!")

	if basis['ptype'] == 'adpt':
		pprop=p_update
	elif basis['ptype'] == 'frzn':
		pprop=frozen
	else:
		sys.exit("Unrecognized keyword in basis ptype!")

	if basis['atype'] == 'adpt':
		aprop=a_update
	elif basis['atype'] == 'frzn':
		aprop=frozen
	else:
		sys.exit("Unrecognized keyword in basis atype!")

	if basis['stype'] == 'adpt':
		print("Adaptable s not implemented yet! Reverting to frozen...")
		sprop=frozen
	elif basis['stype'] == 'frzn':
		sprop=frozen
	else:
		sys.exit("Unrecognized keyword in basis stype!")

	props=[qprop,pprop,aprop,sprop]
	return(props)
