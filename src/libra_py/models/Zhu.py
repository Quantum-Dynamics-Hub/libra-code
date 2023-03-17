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
.. module:: Zhu
   :platform: Unix, Windows
   :synopsis: This module implements various model problems that are used for decoherence testing, as described
              in the work of Chaoyuan Zhu "Restoring electronic coherence/decoherence for a trajectory-based 
              nonadiabatic molecular dynamics" Sci. Rep. 6, 24198 (2016)
.. moduleauthor:: Alexey V. Akimov

"""
import os
import sys
import math
import copy

#from numba import jit

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units



class tmp:
    pass    

def dual_RZD(q, params, full_id):
    """
   
    Dual Rosen-Zener-Demkov model

    Defined in: Chaoyuan Zhu "Restoring electronic coherence/decoherence for a trajectory-based 
                nonadiabatic molecular dynamics" Sci. Rep. 6, 24198 (2016)
    Also see:   J. Xu and L. Wang JCP 150, 164101 (2019) 

    H_00 = 0
    H_11 = A
    H_01 = H_10 = B*[ exp( -C*(x-Z)^2 ) + exp( -C*(x+Z)^2 ) ) ] 
 
    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.025, units: Ha]
            * **params["B"]** ( double ):  [ default: 0.025,  units: Ha]
            * **params["C"]** ( double ):  [ default: 0.7 units: Bohr^-2]
            * **params["Z"]** ( double ):  [ default: 3.0, units: Bohr]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """


    critical_params = [ ] 
    default_params = {"A":0.025, "B":0.025, "C":0.7, "Z":3.0 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    Z = params["Z"]



    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.get(0, indx)


    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    H01, dH01 = 0.0, 0.0

    e_plus = math.exp(-C*(x+Z)*(x+Z))
    e_minus = math.exp(-C*(x-Z)*(x-Z))

    H01 = B * ( e_minus + e_plus)  
    dH01 = -2.0 * B * C * ( (x-Z)*e_minus + (x+Z)*e_plus ) 

    
    Hdia.set(0,0, (0.0+0.0j) );     Hdia.set(0,1, H01*(1.0+0.0j));
    Hdia.set(1,0, H01*(1.0+0.0j));  Hdia.set(1,1, A*(1.0+0.0j));


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




def dual_LZS(q, params, full_id):
    """
   
    Dual Landau-Zener-Stuckelberg model

    Defined in: Chaoyuan Zhu "Restoring electronic coherence/decoherence for a trajectory-based 
                nonadiabatic molecular dynamics" Sci. Rep. 6, 24198 (2016)
    Also see:   J. Xu and L. Wang JCP 150, 164101 (2019) 

    H_00 = 0
    H_11 = E0 - A * exp(-B * x^2)
    H_01 = H_10 = C * exp( -D * x^2 )
 
    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.1, units: Ha]
            * **params["B"]** ( double ):  [ default: 0.28,  units: Bohr^-2]
            * **params["C"]** ( double ):  [ default: 0.01 units: Ha]
            * **params["D"]** ( double ):  [ default: 0.06, units: Bohr^-2]
            * **params["E0"]** ( double ): [ default: 0.03, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """


    critical_params = [ ] 
    default_params = {"A":0.1, "B":0.28, "C":0.01, "D":0.06, "E0":0.03 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]
    E0 = params["E0"]



    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.get(0, indx)
  
    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    H01, dH01 = 0.0, 0.0


    e_b = math.exp(-B*x*x)
    e_d = math.exp(-D*x*x)

    H11 = E0 - A * e_b
    H01 = C * e_d

    dH11 = 2.0 * A * B * e_b * x
    dH01 = -2.0 * C * D * e_d * x

    
    Hdia.set(0,0, (0.0+0.0j) );     Hdia.set(0,1, H01*(1.0+0.0j));
    Hdia.set(1,0, H01*(1.0+0.0j));  Hdia.set(1,1, H11*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 0.0+0.0j );         d1ham_dia[i].set(0,1, dH01*(1.0+0.0j));
        d1ham_dia[i].set(1,0, dH01*(1.0+0.0j));   d1ham_dia[i].set(1,1, dH11*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def Renner_Teller(q, params, full_id):
    """
   
    Renner-Teller

    H_00 = A*(1.0 - exp(-B*x)) x>0,  
         =-A*(1.0 - exp(B*x) ) x<0
    H_11 = -H_00
    H_01 = C * x^2 * exp(-D*x^2)

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.005, units: Ha]
            * **params["B"]** ( double ):  [ default: 1.600, units: Bohr^-1]
            * **params["C"]** ( double ):  [ default: 0.010, units: Ha]
            * **params["D"]** ( double ):  [ default: 1.000, units: Bohr^-2]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """


    critical_params = [ ] 
    default_params = {"A":0.010, "B":1.600, "C":0.005, "D":1.000 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]



    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x = q.get(0, indx)

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    V11,dV11 = 0.0, 0.0
    if x>0:
        e = math.exp(-B*x)
        V11 = A*(1.0 - e)
        dV11 = A*B*e  
    else:
        e = math.exp(B*x)
        V11 = -A*(1.0 - e)
        dV11 = A*B*e

    x2 = x*x
    e = math.exp(-D*x2)
    V = C * x2 * e
    dV = -2.0*x*V + 2.0*x*C*e 

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


