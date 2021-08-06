"""
..module:: QTAG_basis
  :platform: Unix, Windows
  :synopsis: This module contains update functions for the {q,p,a,s} basis components.

..moduleauthors :: Matthew Dutra
"""

from liblibra_core import *
from QTAG_config import univ

"""
.. py:function:: frozen(param_in[, *args])
   Returns the parameter *param_in* unchanged; the 
   non-update function.
"""
 
def frozen(param_in,*args):
	"""The non-update function.
        Args:
                param_in (float): Basis parameter q, p, a, or s.

	Returns:
		param_in (float): The same input parameter q, p, a, or s.
        """

	return(param_in)

"""
.. py:function:: q_adapt(q_in,mom)
   Returns the updated position *q_out* from an input
   position *q_in* and the momentum *mom*, calculated
   via q_out = q_in+mom*dt/m (i.e. symplectic).
"""

def q_adapt(q_in,mom):
	"""Returns the updated position *q_out* from an input position *q_in* and the momentum *mom*, calculated via q_out = q_in+mom*dt/m (i.e. symplectic).

        Args:
                q_in (float): The basis position q to be updated.

		mom (float): Momentum p of basis with position q.
        Returns:
                q_out (float): Updated basis position q.
        """

	q_out = q_in+mom*univ['dt']/univ['mass']
	return(q_out)

"""
.. py:function:: p_adapt(mom)
   Returns the updated basis parameter *p_out* from the 
   momentum *mom*. Currently these two things are equal,
   but they don't necessarily need to be.
"""

def p_adapt(mom):
	"""Returns the updated basis parameter *p_out* from the momentum *mom*. Currently these two things are equal, but they don't necessarily need to be.

        Args:
                mom (float): Basis momentum.

        Returns:
                p_out (float): Updated basis parameter p.
        """

	p_out = 1.0*mom
	return(p_out)

"""
.. py:function:: a_adapt(a_in,gmom)
   Returns the updated basis width *a_out*, calculated
   from the old basis width *a_in* and the momentum
   gradient *gmom*.
"""

def a_adapt(a_in,gmom):
	"""Returns the updated basis width *a_out*, calculated from the old basis width *a_in* and the momentum gradient *gmom*.

        Args:
                a_in (float): Basis width a to be updated.

		gmom (float): Momentum gradient dp/dx at the basis position with width a.
        Returns:
                a_out (float): Updated basis width a.
        """

	an1 = MATRIX(univ['ntraj'],1)
	an1.dot_product(a_in,gmom)
	a_out = a_in-2.0*an1*univ['dt']/univ['mass']
	return(a_out)
	
