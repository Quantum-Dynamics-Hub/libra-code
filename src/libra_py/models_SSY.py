#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file models_SSY.py 
#
# This module implements the 2D, 2-level model of Shenvi-Subotnik-Yang
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



def SSY(q, params):
    """
    Shenvi-Subotnik-Yang, 2-level, 2-dim. problem

    \param[in] q [2 x 1, MATRIX] coordinates of the particles 
    \param[in] params [dictionary] parameters of the model. 

    Should contain the key - value pairs: 
        key          value        description

      "E0" -    (double)      
      "A"  -    (double)
      "B"  -    (double)
      "C"  -    (double)
      "D"  -    (double)

    Ref: Shenvi, N.; Subotnik, J.; Yang, W. JCP 2011, 135, 024101
    
    """

    ndof = q.num_of_rows  # the number of nuclear DOFs

    if ndof != 2:
        print "Error: The SSY Hamiltonian takes 2 nuclear DOFs only. Given = ", ndof
        print "Exiting..."
        sys.exit(0)

    E0 = params["E0"]
    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]


    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in xrange(2):
        obj.d1ham_dia.append( CMATRIX(2,2) )
        obj.dc1_dia.append( CMATRIX(2,2) )


    #=========== Energies & Derivatives ===============
    x, y = q.get(0), q.get(1)

    # H_00
    obj.ham_dia.set(0,0, E0*(-1.0+0.0j))
    obj.d1ham_dia[0].set(0,0, 0.0+0.0j)
    obj.d1ham_dia[1].set(0,0, 0.0+0.0j)

    # H_11
    z = -A*math.exp(-B * (0.75*(x+y)**2 + 0.25*(x-y)**2) )
    dzdx = -B * (1.5*(x+y) + 0.5*(x-y)) * z
    dzdy = -B * (1.5*(x+y) - 0.5*(x-y)) * z
    obj.ham_dia.set(1,1, z * (1.0+0.0j))
    obj.d1ham_dia[0].set(1,1, dzdx * (1.0+0.0j))
    obj.d1ham_dia[1].set(1,1, dzdy * (1.0+0.0j))


    # H_01 and H_10
    z = C*math.exp(-D * (0.25*(x+y)**2 + 0.75*(x-y)**2) )
    dzdx = -D * (0.5*(x+y) + 1.5*(x-y)) * z
    dzdy = -D * (0.5*(x+y) - 1.5*(x-y)) * z

    obj.ham_dia.set(0,1, z * (1.0+0.0j) )
    obj.d1ham_dia[0].set(0,1, dzdx * (1.0+0.0j) )
    obj.d1ham_dia[1].set(0,1, dzdy * (1.0+0.0j) )

    obj.ham_dia.set(1,0, z * (1.0+0.0j) )
    obj.d1ham_dia[0].set(1,0, dzdx * (1.0+0.0j) )
    obj.d1ham_dia[1].set(1,0, dzdy * (1.0+0.0j) )


    return obj



def get_SSY():
    """
    Parameters from:  Shenvi, N.; Subotnik, J.; Yang, W. JCP 2011, 135, 024101

    """

    params = {}
    params["E0"] = 0.05
    params["A"] = 0.15
    params["B"] = 0.14
    params["C"] = 0.015
    params["D"] = 0.06

    return params
   
