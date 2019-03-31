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
.. module:: models_Tully
   :platform: Unix, Windows
   :synopsis: This module implements 3 Tully models for testing NA-MD dynamics
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
import common_utils as comn
import units


class tmp:
    pass    

def Tully1_py(q, params):
    """
   
    Pure Python implementation of the Tully model I = Simple Avoided Crossing (SAC):

    H_00 = A*(1.0-exp(-B*x)) x>0,  
         = A*(exp(B*x)-1.0 ) x<0
    H_11 = -H_00
    H_01 = C*exp(-D*x^2)

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.010, units: Ha]
            * **params["B"]** ( double ):  [ default: 1.600, units: Bohr^-1]
            * **params["C"]** ( double ):  [ default: 0.005, units: Ha]
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
  
    x = q.get(0)
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



def Tully1(q, params):
    """
   
    The implementation that calls the C++ implementation of Tully model I = Simple Avoided Crossing (SAC):

    H_00 = A*(1.0-exp(-B*x)) x>0,  
         = A*(exp(B*x)-1.0 ) x<0
    H_11 = -H_00
    H_01 = C*exp(-D*x^2)

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.010, units: Ha]
            * **params["B"]** ( double ):  [ default: 1.600, units: Bohr^-1]
            * **params["C"]** ( double ):  [ default: 0.005, units: Ha]
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


    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2)
    obj.d1ham_dia = CMATRIXList();  obj.d1ham_dia.append( CMATRIX(2,2) )
    obj.dc1_dia = CMATRIXList();  obj.dc1_dia.append( CMATRIX(2,2) )

    # Convert MATRIX to doubleList()
    qq = doubleList();  qq.append(q.get(0))
    prms = doubleList() # will use the default values
    prms.append(A); prms.append(B); prms.append(C); prms.append(D);
    
    model_SAC(obj.ham_dia, obj.ovlp_dia, obj.d1ham_dia, obj.dc1_dia, qq, prms); 

    return obj



def Tully2(q, params):
    """
   
    The implementation that calls the C++ implementation of Tully model II = Double Avoided Crossing (DAC):

    H_00 = 0.0
    H_11 = E - A*exp(-B*x^2)
    H_01 = C*exp(-D*x^2)

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.100, units: Ha]
            * **params["B"]** ( double ):  [ default: 0.028, units: Bohr^-2]
            * **params["C"]** ( double ):  [ default: 0.015, units: Ha]
            * **params["D"]** ( double ):  [ default: 0.060, units: Bohr^-2]
            * **params["E"]** ( double ):  [ default: 0.050, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [ ] 
    default_params = {"A":0.10, "B":0.28, "C":0.015, "D":0.060, "E":0.050 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]
    E = params["E"]



    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2)
    obj.d1ham_dia = CMATRIXList();  obj.d1ham_dia.append( CMATRIX(2,2) )
    obj.dc1_dia = CMATRIXList();  obj.dc1_dia.append( CMATRIX(2,2) )

    # Convert MATRIX to doubleList()
    qq = doubleList();  qq.append(q.get(0))
    prms = doubleList() # will use the default values
    prms.append(A); prms.append(B); prms.append(C); prms.append(D);  prms.append(E);
    
    model_DAC(obj.ham_dia, obj.ovlp_dia, obj.d1ham_dia, obj.dc1_dia, qq, prms); 

    return obj



def Tully3(q, params):
    """
   
    The implementation that calls the C++ implementation of Tully model III =
    Extended Coupling With Reflection  (ECWR):

    H_00 = A
    H_11 = -H_00
    H_01 = B*exp(C*x);          x <= 0
           B*(2.0 - exp(-C*x)); x > 0

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A"]** ( double ):  [ default: 0.0006, units: Ha ]
            * **params["B"]** ( double ):  [ default: 0.1000, units: Ha ]
            * **params["C"]** ( double ):  [ default: 0.9000, units: Bohr^-1 ]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
 
    """

    critical_params = [ ] 
    default_params = {"A":0.0006, "B":0.1000, "C":0.9000 }
    comn.check_input(params, default_params, critical_params)

    A = params["A"]
    B = params["B"]
    C = params["C"]


    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2)
    obj.d1ham_dia = CMATRIXList();  obj.d1ham_dia.append( CMATRIX(2,2) )
    obj.dc1_dia = CMATRIXList();  obj.dc1_dia.append( CMATRIX(2,2) )

    # Convert MATRIX to doubleList()
    qq = doubleList();  qq.append(q.get(0))
    prms = doubleList() # will use the default values
    prms.append(A); prms.append(B); prms.append(C)
    
    model_ECWR(obj.ham_dia, obj.ovlp_dia, obj.d1ham_dia, obj.dc1_dia, qq, prms); 

    return obj



