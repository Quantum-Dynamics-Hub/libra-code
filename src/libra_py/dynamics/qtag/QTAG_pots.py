"""
..module:: QTAG_pots
  :platform: Unix, Windows
  :synopsis: This module contains functions for computing potential values for various models.

..moduleauthors :: Matthew Dutra
"""

import sys
from liblibra_core import *

#External module for NaFH system, not needed for general QTAG.
#import nafh2dpy

import numpy as np
from QTAG_config import univ, model_params

def model_ho_diab_nD(q,nsurf):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for the coupled harmonic oscillator model in the diabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
                q (MATRIX): The coordinates at which to calculate the potential.

                nsurf (integer): Integer specifying the energetic surface: 1=ground state; 2=excited state; 3=coupling potential

        Returns:
                fx (float): The value of the potential function at the point x.

		dfx (list of floats): The values of the first derivatives of the potential function at the point x, [df/dx1,df/dx2...].

                d2fx (list of floats): The values of the second derivatives of the potential function at the point x, [d2f/d2x1,d2f,d2x2,...].
        """

	ndof=univ['ndof']
	y0=[]
	x0=model_params["x0"];y0.append(5.0*np.sqrt(model_params['k1'][0]))
	for i in range(ndof-1):
		y0.append(0.0)

	fx=0.0;dfx=[];d2fx=[]
	if (nsurf==1):
		for i in range(ndof):
			x=q.get(i)
			fx+=0.5*model_params['k1'][i]*x**2
			dfx.append(model_params['k1'][i]*x)
			d2fx.append(model_params['k1'][i])
	elif (nsurf==2):
		for i in range(ndof):
			x=q.get(i)
			fx+=0.5*model_params['k2'][i]*(x-x0[i])**2+y0[i]
			dfx.append(model_params['k2'][i]*(x-x0[i]))
			d2fx.append(model_params['k2'][i])
	elif (nsurf==3):
		for i in range(ndof):
			x=q.get(i)

			fx+=model_params['d1'][i]*np.exp(-model_params['d2'][i]*(x-model_params['d3'][i])**2)

			dfx.append(-2.0*model_params['d1'][i]*model_params['d2'][i]*(x-model_params['d3'][i])* \
			np.exp(-model_params['d2'][i]*(x-model_params['d3'][i])**2))

			d2fx.append(-2.0*model_params['d1'][i]*model_params['d2'][i]*np.exp(-model_params['d2'][i]* \
			(x-model_params['d3'][i])**2)-dfx[i]*2.0*model_params['d2'][i]*(x-model_params['d3'][i]))

	return(fx,dfx,d2fx)

def model_t1_diab_nD(q,nsurf):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for Tully model I in the diabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
                q (MATRIX): The ndof-by-1 MATRIX of coordinate values at which to calculate the potential.

                nsurf (integer): Integer specifying the energetic surface: 1=ground state; 2=excited state; 3=coupling potential

        Returns:
                fx (float): The value of the potential function at the point x.

                dfx (float): The value of the first derivative of the potential function at the point x.

                d2fx (float): The value of the second derivative of the potential function at the point x.
        """

#	a=0.01;b=1.147
#	d1=0.005;d2=1.0;d3=0.0
	x=q.get(0)

	fx=0.0;dfx=[];d2fx=[]
	if (nsurf==1):
		fx=model_params['a']*(1.0+np.tanh(model_params['b']*x))
		dfx.append(model_params['a']*model_params['b']*(1.0-np.tanh(model_params['b']*x)**2))
		d2fx.append(-2.0*model_params['a']*model_params['b']**2* \
		np.tanh(model_params['b']*x)*(1.0-np.tanh(model_params['b']*x)**2))
	elif (nsurf==2):
		fx=model_params['a']*(1.0-np.tanh(model_params['b']*x))
		dfx.append(-model_params['a']*model_params['b']*(1.0-np.tanh(model_params['b']*x)**2))
		d2fx.append(2.0*model_params['a']*model_params['b']**2* \
		np.tanh(model_params['b']*x)*(1.0-np.tanh(model_params['b']*x)**2))
	elif (nsurf==3):
		fx=model_params['d1']*np.exp(-model_params['d2']*(x-model_params['d3'])**2)

		dfx.append(-2.0*model_params['d1']*model_params['d2']*(x-model_params['d3'])* \
		np.exp(-model_params['d2']*(x-model_params['d3'])**2))

		d2fx.append(-2.0*model_params['d1']*model_params['d2']*np.exp(-model_params['d2']* \
		(x-model_params['d3'])**2)-dfx*2.0*model_params['d2']*(x-model_params['d3']))

	return(fx,dfx,d2fx)

def model_t2_diab_nD(q,nsurf):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for Tully model II in the diabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
                q (MATRIX): The ndof-by-1 MATRIX of coordinate values at which to calculate the potential.

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

def model_t3_adiab_nD(q,nsurf):
	"""Returns the values for the potential (*fx*) and its first (*dfx*) and second (*d2fx*) derivatives at a given point *x* for Tully model III in the adiabatic representation. The parameter *nsurf* specifies which surface is being calculated, and *model_params* is a dictionary containing all relevant parameter values for the model.

        Args:
		q (MATRIX): The ndof-by-1 MATRIX of coordinate values at which to calculate the potential.

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
NaFH test system, but relies on an external module.
def model_nafh_adiab_nD(q,nsurf):
	r=np.zeros([1,3],dtype='d')
	e=np.zeros([1],dtype='d')
	de=np.zeros([3,1],dtype='d')

	if (univ['ndof'] == 1):
		rFH=q.get(0);rNaF=3.779
	elif (univ['ndof'] == 2):
		rFH=q.get(0);rNaF=q.get(1)

	r[0,0]=rNaF+rFH
	r[0,1]=rFH
	r[0,2]=rNaF

	if nsurf == 2:
		nsurf_conv=1
	elif nsurf == 3:
		nsurf_conv=2
	else:
		nsurf_conv=3

	nafh2dpy.prepot()
	nafh2dpy.pot(r,e,de,nsurf_conv,1)

	energy=e[0]
	if (univ['ndof'] == 1):
		denergy=[de[1,0]]
	elif (univ['ndof'] == 2):
		denergy=[de[1,0],de[2,0]]

	return(energy,denergy,0)
"""

def exact_gauss_cpl(qpasi,qpasj,*args,**model_params):
	"""Returns the value *v*=<gi|Vcpl|gj> for a Gaussian coupling potential computed in exact fashion between basis functions on different surfaces defined by the parameters of *qpasi* and *qpasj*.

        Args:
                qpasi (list): List containing the 1-by-4 basis parameter matrices for the basis function at point i.

                qpasj (list): List containing the 1-by-4 basis parameter matrices for the basis function at point j.

        Returns:
                v (complex): The exact value of the Gaussian coupling potential between basis functions i and j.
        """

	v=complex(1.0,0.0)
	q1,p1,a1=qpasi[0], qpasi[1], qpasi[2]
	q2,p2,a2=qpasj[0], qpasj[1], qpasj[2]
	ndof=univ['ndof']

	for i in range(ndof):
		et1=-(q1.get(i)-model_params['d3'][i])**2*a1.get(i)**2/2.0-(q2.get(i)-model_params['d3'][i])**2*a2.get(i)**2/2.0+(p1.get(i)-p2.get(i))**2/2.0
		et2=(q1.get(i)-model_params['d3'][i])*((model_params['d3'][i]-q2.get(i))*a2.get(i)+1j*(p1.get(i)-p2.get(i)))*a1.get(i)
		et3=1j*(q2.get(i)-model_params['d3'][i])*(p1.get(i)-p2.get(i))*a2.get(i)
		et4=(a1.get(i)+2*model_params['d2'][i]+a2.get(i))*(a1.get(i)+a2.get(i))

		v*=np.exp(2.0*model_params['d2'][i]*(et1+et2+et3)/et4)*model_params['d1'][i]* \
		np.sqrt(2.0*a1.get(i)+2.0*a2.get(i))/np.sqrt(2.0*a1.get(i)+4.0*model_params['d2'][i]+2.0*a2.get(i))

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
