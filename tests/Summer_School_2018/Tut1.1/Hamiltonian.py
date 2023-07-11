#*********************************************************************************
#* Copyright (C) 2017-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/


import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


class tmp:
    pass    


def compute_model(q, params, full_id):
    """
    Hdia = 0.5*k*x^2   
    Sdia =  1.0
    Ddia  = 0.0
    """
    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.col(indx).get(0)
    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    obj = tmp()
    obj.ham_dia = CMATRIX(1,1)
    obj.ovlp_dia = CMATRIX(1,1)
    obj.d1ham_dia = CMATRIXList();  obj.d1ham_dia.append( CMATRIX(1,1) )
    obj.dc1_dia = CMATRIXList();  obj.dc1_dia.append( CMATRIX(1,1) )
 
    obj.ham_dia.set(0,0, 0.5*k*x*x*(1.0+0.0j) )
    obj.ovlp_dia.set(0,0, 1.0+0.0j)

    for i in [0]:  # do this for all DOFs
        obj.d1ham_dia[i].set(0,0, k*x*(1.0+0.0j) )  # Put the dU/dx here
        obj.dc1_dia[i].set(0,0, 0.0+0.0j)

    return obj
    
