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
    R = q.col(indx)  # 3*N - dimensional vector
    ndof = R.num_of_rows
    nat = ndof/3

    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    obj = tmp()
    obj.ham_dia = CMATRIX(1,1)
    obj.ovlp_dia = CMATRIX(1,1)
    obj.d1ham_dia = CMATRIXList();  
    obj.dc1_dia = CMATRIXList();  

    for i in xrange(ndof):  # do this for all DOFs
        obj.d1ham_dia.append( CMATRIX(1,1) )
        obj.dc1_dia.append( CMATRIX(1,1) )

     
    Etot = 0.0
    for i in xrange(nat-1):  # do this for all DOFs
        ri = VECTOR(R.get(3*i), R.get(3*i+1), R.get(3*i+2))
        rj = VECTOR(R.get(3*(i+1)), R.get(3*(i+1)+1), R.get(3*(i+1)+2))  # i+1
        rij = ri - rj
        d = rij.length()
        fij = k * (d-x0) * ( rij / d )

        Etot = Etot + 0.5*k*(d-x0)**2        

        obj.d1ham_dia[3*i].add(0,0, fij.x*(1.0+0.0j) )  
        obj.d1ham_dia[3*i+1].add(0,0, fij.y*(1.0+0.0j) )  
        obj.d1ham_dia[3*i+2].add(0,0, fij.z*(1.0+0.0j) )  

        obj.d1ham_dia[3*(i+1)].add(0,0, fij.x*(-1.0+0.0j) )  
        obj.d1ham_dia[3*(i+1)+1].add(0,0, fij.y*(-1.0+0.0j) )  
        obj.d1ham_dia[3*(i+1)+2].add(0,0, fij.z*(-1.0+0.0j) )  


        obj.dc1_dia[3*i].add(0,0, 0.0+0.0j)
        obj.dc1_dia[3*i+1].add(0,0, 0.0+0.0j)
        obj.dc1_dia[3*i+2].add(0,0, 0.0+0.0j)


    obj.ham_dia.set(0,0, Etot*(1.0+0.0j) )
    obj.ovlp_dia.set(0,0, 1.0+0.0j)


    return obj
    
