#*********************************************************************************                     
#* Copyright (C) 2018-2019 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: models_Henon_Heiles
   :platform: Unix, Windows
   :synopsis: This module implements the Henon-Heiles model potential
.. moduleauthor:: Alexey V. Akimov

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

class tmp:
    pass    


def Henon_Heiles(q, params, full_id=None):
    """

    Implementation of the Henon-Heiles potential

    References:
    Sim, E.; Makri, N. Time-Dependent Discrete Variable Representations for Quantum Wave Packet Propagation
    J. Chem. Phys. 1995, 102, 5616-5625

    Args:
        q ( MATRIX(2,1) ): coordinates of the particle, ndof = 2 
        params ( dictionary ): model parameters

            * **params["lam"]** ( double ): lambda parameter [ default: 0.2, units: Ha/Bohr^3]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(1,1) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(1,1) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 2 CMATRIX(1,1) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(1,1) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = [  ] 
    default_params = { "lam":0.2,  }
    comn.check_input(params, default_params, critical_params)

    lam = params["lam"]

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

    x = q.get(0, indx)
    y = q.get(1, indx)

    x2 = x*x
    y2 = y*y
    lam2 = lam*lam

    Hdia.set(0,0,( 0.5*(x2 + y2) + lam*(x*y2 - x*x2/3.0) + lam2*((x2 + y2)**2 / 16.0 ) )*(1.0+0.0j))
    Sdia.set(0,0,1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, ( x + lam*(y2 - x2) + 0.25 * lam2*(x2 + y2)*x )*(1.0+0.0j))
    d1ham_dia[1].set(0,0, ( y + lam*2.0*y*x   + 0.25 * lam2*(x2 + y2)*y )*(1.0+0.0j))

    #  <dia| d/dR_0| dia >
    dc1_dia[0].set(0,0, 0.0+0.0j)
    dc1_dia[1].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


