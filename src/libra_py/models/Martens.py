#*********************************************************************************                     
#* Copyright (C) 2018-2019 Brendan A. Smith, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: models_Martens
   :platform: Unix, Windows
   :synopsis: This module implements Eckart barrier potentials as used by Martens
.. moduleauthor:: Brendan A. Smith, Alexey V. Akimov

"""

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units



class tmp:
    pass    

def model1(q, params, full_id=None):
    """

    The first of the two model potential described by L. Wang, C.C. Martens, 
    and Y. Zheng, "Entangled trajectory molecular dynamics in 
    multidimensional systems: Two-dimensional quantum tunneling through 
    the Eckart barrier" J. Chem. Phys. 137, 34113 (2012)

    This potential has seperable dof

    Args:
        q ( MATRIX(2,1) ): coordinates of the particle, ndof = 2 
        params ( dictionary ): model parameters

            * **params["Va"]** ( double ): barrier height [ default: 0.00625, units: Ha]
            * **params["Vb"]** ( double ): harmonic potential term [ default: 0.0106, units: Ha/Bohr^2]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(1,1) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(1,1) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 2 CMATRIX(1,1) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(1,1) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    # Define potential specific constants
    critical_params = [ ] 
    default_params = {"Va":0.00625, "Vb":0.0106 }
    comn.check_input(params, default_params, critical_params)


    Va = params["Va"]
    Vb = params["Vb"]


    # Hdia and Sdia are ndia x ndia in dimension
    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)

    # d1ham and dc1_dia are ndia x ndia in dimension, but we have nnucl of them
    d1ham_dia = CMATRIXList();
    dc1_dia = CMATRIXList();

    for i in range(0,2):
        d1ham_dia.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) )


    indx = 0
    if full_id !=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]

    x = q.col(indx).get(0)
    y = q.col(indx).get(1)

    x2 = x*x
    y2 = y*y

    # z = sech(2x)
    # z2 = sech^2(2x)    
    z = (2.0 * math.cosh(2.0*x))/(1.0 + math.cosh(4.0*x))
    z2 = z*z

    Hdia.set(0,0,( Va*z2 + 0.5*Vb*y*y )*(1.0+0.0j))
    Sdia.set(0,0,1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, ( ( -4.0*Va*math.tanh(2.0*x)*z2 ) )*(1.0+0.0j))
    d1ham_dia[1].set(0,0, ( Vb*y )*(1.0+0.0j))

    #  <dia| d/dR_0| dia >
    dc1_dia[0].set(0,0, 0.0+0.0j)
    dc1_dia[1].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model2(q, params, full_id=None):
    """

    The second of the two model potential described by L. Wang, C.C. Martens, 
    and Y. Zheng, "Entangled trajectory molecular dynamics in 
    multidimensional systems: Two-dimensional quantum tunneling through 
    the Eckart barrier" J. Chem. Phys. 137, 34113 (2012)

    Specifically, this potential is Model II 

    Hdia = Va*sech^2*(2q1) + 0.5*Vb*[q2 - Vc(q1^2 - 1)]^2

    Args:
        q ( MATRIX(2,1) ): coordinates of the particle, ndof = 2 
        params ( dictionary ): model parameters

            * **params["Va"]** ( double ): barrier height [ default: 0.00625, units: Ha]
            * **params["Vb"]** ( double ): harmonic potential term [ default: 0.0106, units: Ha/Bohr^2]
            * **params["Vc"]** ( double ): coupling term [ default: 0.4, units: Bohr^2]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(1,1) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(1,1) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 2 CMATRIX(1,1) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(1,1) objects ): derivative coupling in the diabatic basis [ zero ]


    """

    # Define potential specific constants
    critical_params = [ ] 
    default_params = {"Va":0.00625, "Vb":0.0106, "Vc":0.4 }
    comn.check_input(params, default_params, critical_params)

    Va = params["Va"]
    Vb = params["Vb"]
    Vc = params["Vc"]


    # Hdia and Sdia are ndia x ndia in dimension
    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)

    # d1ham and dc1_dia are ndia x ndia in dimension, but we have nnucl of them
    d1ham_dia = CMATRIXList();
    dc1_dia = CMATRIXList();

    for i in range(0,2):
        d1ham_dia.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) )

    indx = 0
    if full_id !=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]

    x = q.col(indx).get(0)
    y = q.col(indx).get(1)

    x2 = x*x
    y2 = y*y

    # z = sech(2x)
    # z2 = sech^2(2x)    
    z = (2.0 * math.cosh(2.0*x))/(1.0 + math.cosh(4.0*x))
    z2 = z*z

    Hdia.set(0,0,( Va*z2 + 0.5*Vb*(y+Vc*(x*x-1.0))**2 )*(1.0+0.0j))
    Sdia.set(0,0,1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, ( ( -4.0*Va*math.tanh(2.0*x)*z2 ) + ( 2.0*Vb*Vc*x*(y + Vc*x*x - Vc) ) )*(1.0+0.0j))
    d1ham_dia[1].set(0,0, ( Vb * ( y + Vc*(x*x-1.0) ) ) * (1.0+0.0j))

    #  <dia| d/dR_0| dia >
    dc1_dia[0].set(0,0, 0.0+0.0j)
    dc1_dia[1].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
