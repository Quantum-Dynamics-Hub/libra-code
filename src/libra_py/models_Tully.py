#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file Tully.py 
#
# This module implements 3 Tully models for testing NA-MD dynamics
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

def Tully1(q, params):
    """
    Implementation of the Tully model I

    \param[in] q [ndof x 1, MATRIX] coordinate of the particle
    \param[in] params [dictionary] parameters of the model
 
    """


    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  

    x = q.get(0)
    A = 0.01
    B = 1.6
    C = 0.005
    D = 1.0

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    V11,dV11 = 0.0, 0.0
    if x>0:
        V11 = A*(1.0 - math.exp(-B*x))
        dV11 =  A*B*math.exp(-B*x)
    else:
        V11 = -A*(1.0 - math.exp(B*x))
        dV11 =  A*B*math.exp(B*x)
    
    V = C * math.exp(-D*x*x)
    dV = -2.0*x*C*D*math.exp(-D*x*x)

    Hdia.set(0,0, V11*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));      Hdia.set(1,1, V11*(-1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, dV11*(1.0+0.0j) );  d1ham_dia[i].set(0,1, dV*(1.0+0.0j));
        d1ham_dia[i].set(1,0, dV*(1.0+0.0j));     d1ham_dia[i].set(1,1, dV11*(-1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


