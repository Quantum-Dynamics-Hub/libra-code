"""
..module:: QTAG_pots
  :platform: Unix, Windows
  :synopsis: This module contains functions for computing potential values for various models.

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
.. py:function:: model_ho_diab(x,nsurf,**model_params)
   Returns the values for the potential (*fx*) and its first (*dfx*) and
   second (*d2fx*) derivatives at a given point *x* for the coupled
   harmonic oscillator model in the diabatic representation. The parameter
   *nsurf* specifies which surface is being calculated, and *model_params*
   is a dictionary containing all relevant parameter values for the model.
"""

def model_ho_diab(x,nsurf,**model_params):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for the coupled harmonic oscillator model in the diabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
                x (float): The position at which to calculate the potential.

                nsurf (integer): Integer specifying the energetic surface: 1=ground state; 2=excited state; 3=coupling potential

        Returns:
                fx (float): The value of the potential function at the point x.

		dfx (float): The value of the first derivative of the potential function at the point x.

                d2fx (float): The value of the second derivative of the potential function at the point x.
        """

	x0=1.0;y0=5.0*np.sqrt(model_params['k1']);

	if (nsurf==1):
		fx=0.5*model_params['k1']*x**2
		dfx=model_params['k1']*x
		d2fx=model_params['k1']
	elif (nsurf==2):
		fx=0.5*model_params['k2']*(x-x0)**2+y0
		dfx=model_params['k2']*(x-x0)
		d2fx=model_params['k2']
	elif (nsurf==3):
		fx=model_params['d1']*np.exp(-model_params['d2']*(x-model_params['d3'])**2)

		dfx=-2.0*model_params['d1']*model_params['d2']*(x-model_params['d3'])* \
		np.exp(-model_params['d2']*(x-model_params['d3'])**2)

		d2fx=-2.0*model_params['d1']*model_params['d2']*np.exp(-model_params['d2']* \
		(x-model_params['d3'])**2)-dfx*2.0*model_params['d2']*(x-model_params['d3'])

	return(fx,dfx,d2fx)

"""
.. py:function:: model_t1_diab(x,nsurf,**model_params)
   Returns the values for the potential (*fx*) and its first (*dfx*) and
   second (*d2fx*) derivatives at a given point *x* for Tully model I
   in the diabatic representation. The parameter *nsurf* specifies which 
   surface is being calculated, and *model_params* is a dictionary 
   containing all relevant parameter values for the model.
"""

def model_t1_diab(x,nsurf,**model_params):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for Tully model I in the diabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
                x (float): The position at which to calculate the potential.

                nsurf (integer): Integer specifying the energetic surface: 1=ground state; 2=excited state; 3=coupling potential

        Returns:
                fx (float): The value of the potential function at the point x.

                dfx (float): The value of the first derivative of the potential function at the point x.

                d2fx (float): The value of the second derivative of the potential function at the point x.
        """

#	a=0.01;b=1.147
#	d1=0.005;d2=1.0;d3=0.0
	if (nsurf==1):
		fx=model_params['a']*(1.0+np.tanh(model_params['b']*x))
		dfx=model_params['a']*model_params['b']*(1.0-np.tanh(model_params['b']*x)**2)
		d2fx=-2.0*model_params['a']*model_params['b']**2* \
		np.tanh(model_params['b']*x)*(1.0-np.tanh(model_params['b']*x)**2)
	elif (nsurf==2):
		fx=model_params['a']*(1.0-np.tanh(model_params['b']*x))
		dfx=-model_params['a']*model_params['b']*(1.0-np.tanh(model_params['b']*x)**2)
		d2fx=2.0*model_params['a']*model_params['b']**2* \
		np.tanh(model_params['b']*x)*(1.0-np.tanh(model_params['b']*x)**2)
	elif (nsurf==3):
		fx=model_params['d1']*np.exp(-model_params['d2']*(x-model_params['d3'])**2)

		dfx=-2.0*model_params['d1']*model_params['d2']*(x-model_params['d3'])* \
		np.exp(-model_params['d2']*(x-model_params['d3'])**2)

		d2fx=-2.0*model_params['d1']*model_params['d2']*np.exp(-model_params['d2']* \
		(x-model_params['d3'])**2)-dfx*2.0*model_params['d2']*(x-model_params['d3'])

	return(fx,dfx,d2fx)

"""
.. py:function:: model_t2_diab(x,nsurf,**model_params)
   Returns the values for the potential (*fx*) and its first (*dfx*) and
   second (*d2fx*) derivatives at a given point *x* for Tully model II in 
   the diabatic representation. The parameter *nsurf* specifies which 
   surface is being calculated, and *model_params* is a dictionary 
   containing all relevant parameter values for the model.
"""

def model_t2_diab(x,nsurf,**model_params):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for Tully model II in the diabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
                x (float): The position at which to calculate the potential.

                nsurf (integer): Integer specifying the energetic surface: 1=ground state; 2=excited state; 3=coupling potential

        Returns:
                fx (float): The value of the potential function at the point x.

                dfx (float): The value of the first derivative of the potential function at the point x.

                d2fx (float): The value of the second derivative of the potential function at the point x.
        """

#	a=0.1;b=0.28;e0=0.05
#	d1=0.015;d2=0.06;d3=0.0
	if (nsurf==1):
		fx=0.0
		dfx=0.0
		d2fx=0.0
	elif (nsurf==2):
		fx=-model_params['a']*np.exp(-model_params['b']*x**2)+model_params['e0']
		dfx=2.0*model_params['a']*model_params['b']*x*np.exp(-model_params['b']*x**2)
		d2fx=2.0*model_params['a']*model_params['b']*np.exp(-model_params['b']*x**2)- \
		4.0*model_params['a']*model_params['b']**2*x**2*np.exp(-model_params['b']*x**2)
	elif (nsurf==3):
		fx=model_params['d1']*np.exp(-model_params['d2']*(x-model_params['d3'])**2)

		dfx=-2.0*model_params['d1']*model_params['d2']*(x-model_params['d3'])* \
		np.exp(-model_params['d2']*(x-model_params['d3'])**2)

		d2fx=-2.0*model_params['d1']*model_params['d2']*np.exp(-model_params['d2']* \
		(x-model_params['d3'])**2)-dfx*2.0*model_params['d2']*(x-model_params['d3'])

	return(fx,dfx,d2fx)

"""
.. py:function:: model_t3_adiab(x,nsurf,**model_params)
   Returns the values for the potential (*fx*) and its first (*dfx*) and
   second (*d2fx*) derivatives at a given point *x* for Tully model III in 
   the adiabatic representation. The parameter *nsurf* specifies which 
   surface is being calculated, and *model_params* is a dictionary 
   containing all relevant parameter values for the model.
"""

def model_t3_adiab(x,nsurf,**model_params):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for Tully model III in the adiabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
                x (float): The position at which to calculate the potential.

                nsurf (integer): Integer specifying the energetic surface: 1=ground state; 2=excited state; 3=coupling potential

        Returns:
                fx (float): The value of the potential function at the point x.

                dfx (float): The value of the first derivative of the potential function at the point x.

                d2fx (float): The value of the second derivative of the potential function at the point x.
        """

	if (nsurf!=3):
		if (x<0):
			fx=-np.sqrt(model_params['b']**2*(np.exp(model_params['c']*x))**2+model_params['a']**2)
			dfx=-(model_params['b']**2*(np.exp(model_params['c']*x))*model_params['c']* \
			np.exp(model_params['c']*x))/np.sqrt(model_params['b']**2*(np.exp(model_params['c']*x))**2+ \
			model_params['a']**2)

			t1=(model_params['b']**4*np.exp(model_params['c']*x)**2*(model_params['c']* \
			np.exp(model_params['c']*x))**2)/(model_params['b']**2*np.exp(model_params['c']*x)**2+ \
			model_params['a']**2)**(3/2)
			t2=-model_params['b']**2*(model_params['c']*np.exp(model_params['c']*x))**2/ \
			np.sqrt(model_params['b']**2*np.exp(model_params['c']*x)**2+model_params['a']**2)
			t3=-(model_params['b']**2*np.exp(model_params['c']*x)*model_params['c']**2* \
			np.exp(model_params['c']*x))/np.sqrt(model_params['b']**2*np.exp(model_params['c']*x)**2+ \
			model_params['a']**2)
			d2fx=t1+t2+t3
		else:
			fx=-np.sqrt(model_params['b']**2*(2-np.exp(-model_params['c']*x))**2+model_params['a']**2)
			dfx=-(model_params['b']**2*(2-np.exp(-model_params['c']*x))*model_params['c']* \
			np.exp(-model_params['c']*x))/np.sqrt(model_params['b']**2*(2-np.exp(-model_params['c']*x))**2+ \
			model_params['a']**2)

			t1=(model_params['b']**4*(2-np.exp(-model_params['c']*x))**2*(model_params['c']* \
			np.exp(-model_params['c']*x))**2)/(model_params['b']**2*(2-np.exp(-model_params['c']*x))**2+ \
			model_params['a']**2)**(3/2)
			t2=-model_params['b']**2*(model_params['c']*np.exp(-model_params['c']*x))**2/ \
			np.sqrt(model_params['b']**2*(2-np.exp(-model_params['c']*x))**2+model_params['a']**2)
			t3=-(model_params['b']**2*(2-np.exp(-model_params['c']*x))*(-model_params['c']**2)* \
			np.exp(-model_params['c']*x))/np.sqrt(model_params['b']**2*(2-np.exp(-model_params['c']*x))**2+ \
			model_params['a']**2)
			d2fx=t1+t2+t3

		if (nsurf==1):
			fx=-fx
			dfx=-dfx
			d2fx=-d2fx

	else:
		if (x<0):
			fx=model_params['a']*model_params['b']*model_params['c']*np.exp(model_params['c']*x)/ \
			(2*model_params['b']**2*np.exp(2*model_params['c']*x)+2*model_params['a']**2)

			t1=model_params['b']*model_params['c']**2*np.exp(model_params['c']*x)/(2*model_params['a']* \
			((model_params['b']*np.exp(model_params['c']*x))**2/model_params['a']**2+1))
			t2=-(model_params['b']*model_params['c']*np.exp(model_params['c']*x))**2*(model_params['b']* \
			np.exp(model_params['c']*x))/(model_params['a']**3*((model_params['b']* \
			np.exp(model_params['c']*x))**2/model_params['a']**2+1)**2)
			dfx=t1+t2

			t1=model_params['b']*model_params['c']**3*np.exp(model_params['c']*x)/(2*model_params['a']* \
			((model_params['b']*np.exp(model_params['c']*x))**2/model_params['a']**2+1))
			t2=-3*(model_params['b']*model_params['c']**2*np.exp(model_params['c']*x))*(model_params['b']* \
			np.exp(model_params['c']*x))*(model_params['b']*model_params['c']*np.exp(model_params['c']*x))/ \
			(model_params['a']**3*((model_params['b']*np.exp(model_params['c']*x))**2/model_params['a']**2+1)**2)
			t3=-(model_params['b']*model_params['c']*np.exp(model_params['c']*x))**3/(model_params['a']**3* \
			((model_params['b']*np.exp(model_params['c']*x))**2/model_params['a']**2+1)**2)
			t4=4*(model_params['b']*model_params['c']*np.exp(model_params['c']*x))**3*(model_params['b']* \
			np.exp(model_params['c']*x))**2/(model_params['a']**5*((model_params['b']* \
			np.exp(model_params['c']*x))**2/model_params['a']**2+1)**3)
			d2fx=t1+t2+t3+t4
		else:	
			fx=model_params['a']*model_params['b']*model_params['c']*np.exp(-model_params['c']*x)/ \
			(2*model_params['b']**2*np.exp(-2*model_params['c']*x)-8*model_params['b']**2* \
			np.exp(-model_params['c']*x)+2*model_params['a']**2+8*model_params['b']**2)

			t1=-model_params['b']*model_params['c']**2*np.exp(-model_params['c']*x)/(2*model_params['a']* \
			((model_params['b']*(2-np.exp(-model_params['c']*x)))**2/model_params['a']**2+1))
			t2=-(model_params['b']*model_params['c']*np.exp(-model_params['c']*x))**2*(model_params['b']* \
			(2-np.exp(-model_params['c']*x)))/(model_params['a']**3*((model_params['b']* \
			(2-np.exp(-model_params['c']*x)))**2/model_params['a']**2+1)**2)
			dfx=t1+t2

			t1=model_params['b']*model_params['c']**3*np.exp(-model_params['c']*x)/(2*model_params['a']* \
			((model_params['b']*(2-np.exp(-model_params['c']*x)))**2/model_params['a']**2+1))
			t2=-3*(-model_params['b']*model_params['c']**2*np.exp(-model_params['c']*x))*(model_params['b']* \
			(2-np.exp(-model_params['c']*x)))*(model_params['b']*model_params['c']*np.exp(-model_params['c']*x))/ \
			(model_params['a']**3*((model_params['b']*(2-np.exp(-model_params['c']*x)))**2/model_params['a']**2+1)**2)
			t3=-(model_params['b']*model_params['c']*np.exp(-model_params['c']*x))**3/(model_params['a']**3* \
			((model_params['b']*(2-np.exp(-model_params['c']*x)))**2/model_params['a']**2+1)**2)
			t4=4*(model_params['b']*model_params['c']*np.exp(-model_params['c']*x))**3*(model_params['b']* \
			(2-np.exp(-model_params['c']*x)))**2/(model_params['a']**5*((model_params['b']* \
			(2-np.exp(-model_params['c']*x)))**2/model_params['a']**2+1)**3)
			d2fx=t1+t2+t3+t4

	return(fx,dfx,d2fx)

"""
.. py:function:: exact_gauss_cpl(qpasi,qpasj,*args,**model_params)
   Returns the value *v*=<gi|Vcpl|gj> for a Gaussian coupling potential 
   computed in exact fashion between basis functions on different surfaces
   defined by the parameters of *qpasi* and *qpasj*.
"""

def exact_gauss_cpl(qpasi,qpasj,*args,**model_params):
	"""Returns the value *v*=<gi|Vcpl|gj> for a Gaussian coupling potential computed in exact fashion between basis functions on different surfaces defined by the parameters of *qpasi* and *qpasj*.

        Args:
                qpasi (MATRIX): The 1-by-4 basis parameter matrix for the basis function at point i.

                qpasj (MATRIX): The 1-by-4 basis parameter matrix for the basis function at point j.

        Returns:
                v (complex): The exact value of the Gaussian coupling potential between basis functions i and j.
        """

	v=complex(0.0,0.0)
	q1,p1,a1=qpasi.get(0), qpasi.get(1), qpasi.get(2)
	q2,p2,a2=qpasj.get(0), qpasj.get(1), qpasj.get(2)

	et1=-(q1-model_params['d3'])**2*a1**2/2.0-(q2-model_params['d3'])**2*a2**2/2.0+(p1-p2)**2/2.0
	et2=(q1-model_params['d3'])*((model_params['d3']-q2)*a2+1j*(p1-p2))*a1
	et3=1j*(q2-model_params['d3'])*(p1-p2)*a2
	et4=(a1+2*model_params['d2']+a2)*(a1+a2)

	v=np.exp(2.0*model_params['d2']*(et1+et2+et3)/et4)*model_params['d1']* \
	np.sqrt(2.0*a1+2.0*a2)/np.sqrt(2.0*a1+4.0*model_params['d2']+2.0*a2)

	return(v)

def no_coupling(*args,**model_params):
	"""Returns the value v=<gi|Vcpl|gj>=0 for a system with no coupling potential.

	Args:
		None.

        Returns:
                v (complex): The exact value for no coupling between surfaces (=0+0j).
        """

	v=complex(0.0,0.0)
	return(v)
