#*********************************************************************************                     
#* Copyright (C) 2022 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: Subotnik
   :platform: Unix, Windows
   :synopsis: This module implements various model problems from Joe Subotnik
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
import libra_py.units as units


class tmp:
    pass    

def dumbbell_geometry(q, params, full_id=None):
    """
   
    Dubmbell Geometry = a symmetrized version of the Tully's ECWR model. 

    Defined in: J. E. Subotnik and N. Shenvi JCP 134, 024105 (2011)
    Also see:   J. Xu and L. Wang JCP 150, 164101 (2019) - the definition

    H_00 = A
    H_11 = -A
    H_01 = H_10 = B [ exp( C*(x-Z) ) + (2 - exp( C*(x+Z) ) ) ]        x < -Z
                = B [ exp( C*(x-Z) ) + exp( -C*(x+Z)  ) ]        -Z < x < Z
                = B [ exp( -C*(x+Z) ) + (2 - exp( -C*(x-Z) ) ) ]      x > Z
 
    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.0006, units: Ha]
            * **params["B"]** ( double ):  [ default: 0.100,  units: Ha]
            * **params["C"]** ( double ):  [ default: 0.9 units: Bohr^-1]
            * **params["Z"]** ( double ):  [ default: 10.000, units: Bohr]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """


    critical_params = [ ] 
    default_params = {"A":0.0006, "B":0.1, "C":0.9, "Z":10.000 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    Z = params["Z"]



    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  
    indx = 0
    if full_id !=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]

    x = q.get(0, indx)

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    H01, dH01 = 0.0, 0.0
    if x<-Z:
        e_plus = math.exp(C*(x+Z))
        e_minus = math.exp(C*(x-Z))

        H01 = B * ( e_minus + (2.0 - e_plus) )
        dH01 = B * C * ( e_minus - e_plus) 

    elif x>-Z and x<Z:
        e_plus = math.exp(-C*(x+Z))
        e_minus = math.exp(C*(x-Z))

        H01 = B * ( e_minus + e_plus) 
        dH01 = B * C * ( e_minus - e_plus) 
         
    else:
        e_plus = math.exp(-C*(x+Z))
        e_minus = math.exp(-C*(x-Z))

        H01 = B * ( e_plus + (2.0 - e_minus) )
        dH01 = -B * C * ( e_plus - e_minus) 

    
    Hdia.set(0,0, A*(1.0+0.0j) );   Hdia.set(0,1, H01*(1.0+0.0j));
    Hdia.set(1,0, H01*(1.0+0.0j));  Hdia.set(1,1, A*(-1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 0.0+0.0j );         d1ham_dia[i].set(0,1, dH01*(1.0+0.0j));
        d1ham_dia[i].set(1,0, dH01*(1.0+0.0j));   d1ham_dia[i].set(1,1, 0.0+0.0j);

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def double_arch_geometry(q, params, full_id=None):
    """
   
    Dubmbell Geometry = a symmetrized version of the Tully's ECWR model. 

    Defined in: J. E. Subotnik and N. Shenvi JCP 134, 024105 (2011)
    Also see:   J. Xu and L. Wang JCP 150, 164101 (2019) - the definition

    H_00 = A
    H_11 = -A
    H_01 = H_10 = B [ -exp( C*(x-Z) ) + exp( C*(x+Z) )  ]             x <-Z
                = B [ -exp( C*(x-Z) ) - exp( -C*(x+Z) ) + 2.0 ]  -Z < x < Z
                = B [ exp( -C*(x-Z) ) - exp( -C*(x+Z) )  ]            x > Z
 
    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.0006, units: Ha]
            * **params["B"]** ( double ):  [ default: 0.100,  units: Ha]
            * **params["C"]** ( double ):  [ default: 0.9 units: Bohr^-1]
            * **params["Z"]** ( double ):  [ default: 4.000, units: Bohr]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """


    critical_params = [ ] 
    default_params = {"A":0.0006, "B":0.1, "C":0.9, "Z":4.000 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    Z = params["Z"]



    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )

    indx = 0
    if full_id !=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]
  
    x = q.get(0, indx)
    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

#    H_01 = H_10 = B [ -exp( C*(x-Z) ) + exp( C*(x+Z) )  ]             x <-Z
#                = B [ -exp( C*(x-Z) ) - exp( -C*(x+Z) ) + 2.0 ]  -Z < x < Z
#                = B [ exp( -C*(x-Z) ) - exp( -C*(x+Z) )  ]            x > Z


    H01, dH01 = 0.0, 0.0
    if x<-Z:
        e_plus = math.exp(C*(x+Z))
        e_minus = math.exp(C*(x-Z))

        H01 = B * ( -e_minus + e_plus)
        dH01 = B * C * ( -e_minus + e_plus) 

    elif x>-Z and x<Z:
        e_plus = math.exp(-C*(x+Z))
        e_minus = math.exp(C*(x-Z))

        H01 = B * ( -e_minus - e_plus + 2.0) 
        dH01 = B * C * ( -e_minus + e_plus) 
         
    else:
        e_plus = math.exp(-C*(x+Z))
        e_minus = math.exp(-C*(x-Z))

        H01 = B * ( e_minus - e_plus) 
        dH01 = -B * C * ( e_minus - e_plus) 

    
    Hdia.set(0,0, A*(1.0+0.0j) );   Hdia.set(0,1, H01*(1.0+0.0j));
    Hdia.set(1,0, H01*(1.0+0.0j));  Hdia.set(1,1, A*(-1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 0.0+0.0j );         d1ham_dia[i].set(0,1, dH01*(1.0+0.0j));
        d1ham_dia[i].set(1,0, dH01*(1.0+0.0j));   d1ham_dia[i].set(1,1, 0.0+0.0j);

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



