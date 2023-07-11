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
.. module:: Granucci_Persico
   :platform: Unix, Windows
   :synopsis: This module implements model Hamiltonians described in
              Granucci, G.; Persico, M.; Zoccante, A. "Including quantum decoherence in surface hopping" 
              J. Chem. Phys. 133, 134111, (2010)
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

def model1(q, params):
    """   
    Originally defined in:
    Granucci, G.; Persico, M. "Critical appraisal of the fewest switches algorithm for surface hopping" 
              J. Chem. Phys. 126, 134114, (2007)   

    H_00 = a1 * exp(-alp1 * x)  + dE
    H_11 = a2 * exp(-alp2 * x)
    H_01 = H_10 = b * exp(-beta * (x-x_c)^2 ) + gamma * sin^2(x)
 
    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["a1"]** ( double ):  [ default: 0.06, units: Ha]
            * **params["a2"]** ( double ):  [ default: 0.5, units: Ha]
            * **params["alp1"]** ( double ):  [ default: 0.2, units: Bohr^-1]
            * **params["alp2"]** ( double ):  [ default: 0.5, units: Bohr^-1]
            * **params["dE"]** ( double ):  [ default: 0.03, units: Ha]
            * **params["b"]** ( double ):  [ default: 0.013, units: Ha]
            * **params["beta"]** ( double ):  [ default: 0.2, units: Bohr^-1]
            * **params["gamma"]** ( double ):  [ default: 0.0002, units: Ha]
            * **params["x_c"]** ( double ):  [ default: 5.0, units: Bohr]


    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [ ] 
    default_params = {"a1":0.006, "a2":0.5, "alp1":0.2, "alp2":0.5,
                      "dE":0.03, "b":0.013, "beta":0.2, "gamma":0.0002, "x_c":5.0 }
    comn.check_input(params, default_params, critical_params)

    a1 = params["a1"]
    a2 = params["a2"]
    alp1 = params["alp1"]
    alp2 = params["alp2"]
    dE = params["dE"]
    b = params["b"]
    beta = params["beta"]
    gamma = params["gamma"]
    x_c = params["x_c"]


    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  
    x = q.get(0)
    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    e1 = a1 * math.exp(-alp1 * x)
    H_00 = e1  + dE
    dH_00 = -alp1 * e1

    e2 = a2 * math.exp(-alp2 * x)
    H_11 = e2
    dH_11 = -alp2 * e2

    e3 = b * math.exp(-beta * (x-x_c)**2 )
    si = math.sin(x)
    cs = math.cos(x)
    H_01 = e3 + gamm * si**2
    dH_01 = -2.0 * beta * (x-x_c) * e3 + gamma * 2.0 * si * cs

    
    Hdia.set(0,0, H_00*(1.0+0.0j) );     Hdia.set(0,1, H_01*(1.0+0.0j));
    Hdia.set(1,0, H_01*(1.0+0.0j));      Hdia.set(1,1, H_11*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, dH_00*(1.0+0.0j) );   d1ham_dia[i].set(0,1, dH_01*(1.0+0.0j) );
        d1ham_dia[i].set(1,0, dH_01*(1.0+0.0j) );   d1ham_dia[i].set(1,1, dH_11*(1.0+0.0j) );

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj





def model2(q, params):
    """   
    This is a 2D conical intersection model originally defined in:
    Granucci, G.; Persico, M.; Zoccante, A. "Including quantum decoherence in surface hopping" 
              J. Chem. Phys. 133, 134111, (2010)

    H_00 = D1 * [exp(-2 * alp1 * (x-x1)) - 2 * exp(-alp1 * (x-x1)) ]  + delta1 +  0.5 * K * y**2
    H_11 = D2 * [exp(-2 * alp2 * (x-x2)) - 2 * exp(-alp2 * (x-x2)) ]  + delta2 +  0.5 * K * y**2
    H_01 = H10 = gamma * y * exp(-beta1 * (x-x3)**2 - beta2 * y**2 )
 
    Args: 
        q ( MATRIX(2,1) ): coordinates of the particle, ndof = 2
        params ( dictionary ): model parameters

            * **params["D1"]** ( double ):  [ default: 0.015, units: Ha]
            * **params["D2"]** ( double ):  [ default: 0.11, units: Ha]
            * **params["delta1"]** ( double ):  [ default: 0.05, units: Ha]
            * **params["delta2"]** ( double ):  [ default: 0.11, units: Ha]
            * **params["alp1"]** ( double ):  [ default: 1.0, units: Bohr^-1]
            * **params["alp2"]** ( double ):  [ default: 0.674, units: Bohr^-1]
            * **params["beta1"]** ( double ):  [ default: 0.5, units: Bohr^-2]
            * **params["beta2"]** ( double ):  [ default: 1.5, units: Bohr^-2]
            * **params["gamma"]** ( double ):  [ default: 2.2295e-2, units: Ha]
            * **params["K"]** ( double ):  [ default: 0.09, units: Ha*Bohr^-2]
            * **params["x1"]** ( double ):  [ default: 3.9, units: Bohr]
            * **params["x2"]** ( double ):  [ default: 3.0, units: Bohr]
            * **params["x3"]** ( double ):  [ default: 5.0, units: Bohr]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 2 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """
 
    critical_params = [ ] 
    default_params = {"D1":0.015, "D2":0.11, "delta1":0.05, "delta2":0.11,
                      "alp1":1.0, "alp2":0.674, "beta1":0.5, "beta2":1.5, "gamma":2.2295e-2,
                      "K":0.09, "x1":3.9, "x2":3.0, "x3":5.0 }
    comn.check_input(params, default_params, critical_params)

    D1 = params["D1"]
    D2 = params["D2"]
    alp1 = params["alp1"]
    alp2 = params["alp2"]
    beta1 = params["beta1"]
    beta2 = params["beta2"]
    gamma = params["gamma"]
    delta1 = params["delta1"]
    delta2 = params["delta2"]
    K = params["K"]
    x1 = params["x1"]
    x2 = params["x2"]
    x3 = params["x3"]

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  
    x, y = q.get(0), q.get(1)
    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);


    e1 = math.exp(-alp1*(x-x1))
    e2 = math.exp(-alp2*(x-x2))
    e3 = math.exp( -beta1*(x-x3)**2 - beta2*y**2)

    H_00 = D1 * (e1**2 - 2.0 * e1 ) + delta1 +  0.5 * K * y**2
    H_11 = D2 * (e2**2 - 2.0 * e2 ) + delta2 +  0.5 * K * y**2
    H_01 = gamma * y * e3
    
    Hdia.set(0,0, H_00*(1.0+0.0j) );     Hdia.set(0,1, H_01*(1.0+0.0j));
    Hdia.set(1,0, H_01*(1.0+0.0j));      Hdia.set(1,1, H_11*(1.0+0.0j));

    # Derivatives
    #  d Hdia / dR_0
    de1 = -alp1 * e1
    dH_00 = 2.0 * D1 * (e1 - 1.0) * de1

    de2 = -alp2 * e2
    dH_11 = 2.0 * D2 * (e2 - 1.0) * de2

    dH_01 = -2.0*beta1*(x-x3)* H_01

    d1ham_dia[0].set(0,0, dH_00*(1.0+0.0j) );   d1ham_dia[0].set(0,1, dH_01*(1.0+0.0j) );
    d1ham_dia[0].set(1,0, dH_01*(1.0+0.0j) );   d1ham_dia[0].set(1,1, dH_11*(1.0+0.0j) );

    #  d Hdia / dR_1
    dH_00 = K * y
    dH_11 = dH_00
    dH_01 = gamma * e3 - 2.0*beta*y* H_01

    d1ham_dia[0].set(0,0, dH_00*(1.0+0.0j) );   d1ham_dia[0].set(0,1, dH_01*(1.0+0.0j) );
    d1ham_dia[0].set(1,0, dH_01*(1.0+0.0j) );   d1ham_dia[0].set(1,1, dH_11*(1.0+0.0j) );


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia   # <dia| d/dR_0| dia >

    return obj

