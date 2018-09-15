#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file models_Henon_Heiles.py 
#
# This module implements the Henon-Heiles model potential
#
#
import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

class tmp:
    pass    


def Henon_Heiles(q, params):
    """
    Implementation of the Henon-Heiles potential

    \param[in] q [ndof x 1, MATRIX] coordinate of the particle, ndof = 2 
    \param[in] params [dictionary] parameters of the model

    References:
    Sim, E.; Makri, N. Time-Dependent Discrete Variable Representations for Quantum Wave Packet Propagation. 
    J. Chem. Phys. 1995, 102, 5616–5625.
 
    """

    lam = 0.2   

    # Hdia and Sdia are ndia x ndia in dimension
    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)

    # d1ham and dc1_dia are ndia x ndia in dimension, but we have nnucl of them
    d1ham_dia = CMATRIXList();
    dc1_dia = CMATRIXList();

    for i in xrange(2):
        d1ham_dia.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) )

    x = q.get(0)
    y = q.get(1)

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


