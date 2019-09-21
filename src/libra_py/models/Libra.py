#*********************************************************************************                     
#* Copyright (C) 2018-2019 Brendan A. Smith, Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
"""
.. module:: models_Libra
   :platform: Unix, Windows
   :synopsis: The Libra models - well, they are not necessarily introduced here for the 
       first time (the special cases might have been already used in different context),
       but these are the models we define for the internal test purposes

.. moduleauthor:: Brendan A. Smith, Alexey V. Akimov

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

    Essentially the spin-boson (Marcus) model
    

              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states
                [ default: 1.0, units: Bohr ]
            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha/Bohr^2]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 


    """

    critical_params = [ ] 
    default_params = {"x0":1.0, "k":0.01, "D":0.0, "V":0.005 }
    comn.check_input(params, default_params, critical_params)

    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  

    x = q.get(0)


    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def model1a(Hdia, Sdia, d1ham_dia, dc1_dia, q, params):
    """

    Same as ::funct:```model1``` just different interface
    

              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    Args: 
        Hdia ( CMATRIX(2,2) ): diabatic Hamiltonian - updated by this function
        Sdia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states - updated by this function [ identity ] 
        d1ham_dia ( list of 1 CMATRIX(2,2) objects ): derivatives of the diabatic Hamiltonian w.r.t. 
            the nuclear coordinate - updated by this function
        dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis - updated 
            by this function [ zero ]
        q ( double ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states
                [ default: 1.0, units: Bohr ]
            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha/Bohr^2]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]

    Returns:       
        None
 
    """

    x = q
    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]


    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);



def convert_Landry_Subotnik2model1(m, Er, w, eps):
    """

    This function converts the parameter of the Landry-Subotnik Hamiltonian
    to those of the ```model1``` form. This mapping is up to the arbitratry energy
    scale shift factor. 

    Args:
        m ( double ): mass of the particle [a.u. of mass ]
        Er ( double ): reorganization energy [ Ha ]
        w ( double ): oscillator frequency [ a.u.^-1 ]
        eps ( double ): driving energy [ Ha ]

    Returns:
        tuple: (k, x0, D): where:
 
            * k ( double ): force constant of the oscillator [ units: Ha*Bohr^-2 ]
            * x0 ( double ): displacement of the energy minima positions [ units: Bohr ]
            * D ( double ): the energy shift of the acceptor with respect to donor [ units: Ha]


    Ref: Landry, B. R.; Subotnik, J. E. J. Chem. Phys. 2011, 135, 191101


    H00 = (1/2) * m * w^2 * x^2  + M x =  (1/2) * m * w^2 * ( x^2  + 2Mx/(m*w^2) ) =
    = (1/2) * m * w^2 * [( x  + M/(m*w^2) )^2 - M^2 /(m^2 *w^4)  ] = 
    = (1/2) * m * w^2 * [( x  + M/(m*w^2) )^2 ]  - M^2 /(2*m*w^2)  

    H11 = (1/2) * m * w^2 * x^2  - M x - e =  (1/2) * m * w^2 * ( x^2  - 2Mx/(m*w^2) ) - e =
    = (1/2) * m * w^2 * [( x  - M/(m*w^2) )^2 - M^2 /(m^2 *w^4)  ] - e = 
    = (1/2) * m * w^2 * [( x  - M/(m*w^2) )^2 ] - M^2 /(2*m*w^2)  - e


    M = sqrt(Er * m*w^2 / 2) =>  2* M^2 / (m*w^2) = Er

    So the connection to the model1 parameters:

    k = (1/2) * m * w^2 

    x0 = 2 * M/(m*w^2)  = M/k = sqrt( Er / k )

    D = -e
    
    """

    k = 0.5*m*w*w
    x0 = math.sqrt(Er / k)
    D = -eps

    return k, x0, D


def get_Landry_Subotnik_set1(gamma_i, eps_i):
    """

    Returns one of the data point in the parameters set from 
    Ref: Landry, B. R.; Subotnik, J. E. J. Chem. Phys. 2011, 135, 191101

    The data set follows their Figure 2

    Args:
        gamma_i ( int ): index of the gamma parameter
        eps_i ( int ): index of the epsilon parameter

    Returns:
        dictionary: parameters in the format suitable for use with model1:

            * **params["gamma"]** ( double ): friction coefficient [ units: (a.u. time)* Ha/amu*Bohr^2 ]
            * **params["V"]** ( double ): electronic coupling between these diabats [ units: Ha]
            * **params["mass"]** ( double ): mass of the particle [ units: amu ]
            * **params["kT"]** ( double ): system's temperature [ units: Ha ]
            * **params["k"]** ( double ): force constante [ units: Ha/Bohr]
            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states [ units: Bohr ]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ units: Ha]

    """

    gammas = [1.875e-5, 3.75e-5, 7.5e-5, 1.5e-4, 3.0e-4, 6.0e-4, 1.2e-3, 2.4e-3]
    
    m = 1.0
    Er = 2.39e-2
    w = 3.5e-4
    eps = 0.015 + eps_i*(0.03 - 0.015)/16.0 

    k, x0, D =  convert_Landry_Subotnik2model1(m, Er, w, eps)


    params =  {}
    params["gamma"] = gammas[gamma_i]
    params["V"] = 5.0e-5
    params["mass"] = 1.0
    params["kT"] = 9.5e-4

    params["k"] = k
    params["x0"] = x0
    params["D"] = D
    
    return params


def get_Landry_Subotnik_set2(V_i):
    """

    Returns one of the data point in the parameters set from 
    Ref: Landry, B. R.; Subotnik, J. E. J. Chem. Phys. 2011, 135, 191101

    The data set follows their Figure 3

    Args:
        V_i ( int ): index of the coupling parameter

    Returns:
        dictionary: parameters in the format suitable for use with model1:

            * **params["gamma"]** ( double ): friction coefficient [ units: (a.u. time)* Ha/amu*Bohr^2 ]
            * **params["V"]** ( double ): electronic coupling between these diabats [ units: Ha]
            * **params["mass"]** ( double ): mass of the particle [ units: amu ]
            * **params["kT"]** ( double ): system's temperature [ units: Ha ]
            * **params["k"]** ( double ): force constante [ units: Ha/Bohr]
            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states [ units: Bohr ]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ units: Ha]

    """

        
    m = 1.0
    Er = 2.39e-2
    w = 3.5e-4
    eps = 0.015

    k, x0, D =  convert_Landry_Subotnik2model1(m, Er, w, eps)


    params =  {}
    params["gamma"] = 0.0024
    params["V"] = math.exp( math.log(1.49e-5) + V_i * (math.log(2.28e-4) - math.log(1.49e-5))/16.0 )
    params["mass"] = 1.0
    params["kT"] = 9.5e-4

    params["k"] = k
    params["x0"] = x0
    params["D"] = D
    
    return params




def model2(q, params):
    """

    Essentially the spin-boson (Marcus) model    

              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D


    Sdia =  I

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia = 0.0


    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states
                [ default: 1.0, units: Bohr ]
            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha/Bohr^2]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]
            * **params["NAC"]** ( double ): NAC in the diabatic basis  [ default: -0.1, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [ ] 
    default_params = {"x0":1.0, "k":0.01, "D":0.0, "V":0.005, "NAC":-0.1 }
    comn.check_input(params, default_params, critical_params)
    x0,k,D,V, nac = params["x0"], params["k"], params["D"], params["V"], params["NAC"]


    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);          dc1_dia[i].set(0,1,nac*(1.0+0.0j));
        dc1_dia[i].set(1,0, nac*(-1.0+0.0j));   dc1_dia[i].set(1,1, 0.0+0.0j);


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def model2a(Hdia, Sdia, d1ham_dia, dc1_dia, q, params):
    """

    Same as ::funct:```model2``` just different interface
    

              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia = 0.0

    Args: 
        Hdia ( CMATRIX(2,2) ): diabatic Hamiltonian - updated by this function
        Sdia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states - updated by this function [ identity ] 
        d1ham_dia ( list of 1 CMATRIX(2,2) objects ): derivatives of the diabatic Hamiltonian w.r.t. 
            the nuclear coordinate - updated by this function
        dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis - updated 
            by this function [ zero ]
        q ( double ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states
                [ default: 1.0, units: Bohr ]
            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha/Bohr^2]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]
            * **params["NAC"]** ( double ): NAC in the diabatic basis  [ default: -0.1, units: Ha]

    Returns:       
        None


    """

    critical_params = [ ] 
    default_params = {"x0":1.0, "k":0.01, "D":0.0, "V":0.005, "NAC":-0.1 }
    comn.check_input(params, default_params, critical_params)
    x0,k,D,V, nac = params["x0"], params["k"], params["D"], params["V"], params["NAC"]

    x = q

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);          dc1_dia[i].set(0,1, nac*(1.0+0.0j));
        dc1_dia[i].set(1,0, nac*(-1.0+0.0j));   dc1_dia[i].set(1,1, 0.0+0.0j);



def model3(q, params):
    """

    More complex spin-boson model    

              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D


    Sdia =    1                          B*exp(-(x-0.5*x0)^2)
             B*exp(-(x-0.5*x0)^2)               1

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia !=0.0


    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states
                [ default: 1.0, units: Bohr ]
            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha/Bohr^2]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]
            * **params["B"]** ( double ): parameter controlling the overlap of the diabatic states
                [ default: 0.05, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]


    """
    critical_params = [ ] 
    default_params = {"x0":1.0, "k":0.01, "D":0.0, "V":0.005, "B":0.05 }
    comn.check_input(params, default_params, critical_params)
    x0,k,D,V,B = params["x0"], params["k"], params["D"], params["V"], params["B"]


    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  

    x = q.get(0)

    ex = B*math.exp(-(x-0.5*x0)**2)*(1.0+0.0j)

    Sdia.set(0,0, 1.0+0.0j);     Sdia.set(0,1, ex );
    Sdia.set(1,0, ex);           Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );  d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);             d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >  = 0.5 * dS/dR
        d = -(x-0.5*x0)*ex
        dc1_dia[i].set(0,0, 0.0+0.0j);  dc1_dia[i].set(0,1, d);
        dc1_dia[i].set(1,0, d);         dc1_dia[i].set(1,1, 0.0+0.0j);


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def model3a(Hdia, Sdia, d1ham_dia, dc1_dia, q, params):
    """

    Same as ::funct:```model3``` just different interface
    

              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =    1                           B*exp(-(x-0.5*x0)^2)
             B*exp(-(x-0.5*x0)^2)               1

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia !=0.0

    Args: 
        Hdia ( CMATRIX(2,2) ): diabatic Hamiltonian - updated by this function
        Sdia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states - updated by this function [ identity ] 
        d1ham_dia ( list of 1 CMATRIX(2,2) objects ): derivatives of the diabatic Hamiltonian w.r.t. 
            the nuclear coordinate - updated by this function
        dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis - updated 
            by this function [ zero ]
        q ( double ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["x0"]** ( double ): displacement of the minimum of one of the diabatic states
                [ default: 1.0, units: Bohr ]
            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha/Bohr^2]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]
            * **params["B"]** ( double ): parameter controlling the overlap of the diabatic states
                [ default: 0.05, units: Ha]

    Returns:       
        None

    """
    critical_params = [ ] 
    default_params = {"x0":1.0, "k":0.01, "D":0.0, "V":0.005, "B":0.05 }
    comn.check_input(params, default_params, critical_params)
    x0,k,D,V,B = params["x0"], params["k"], params["D"], params["V"], params["B"]

    x = q

    ex = B*math.exp(-(x-0.5*x0)**2)*(1.0+0.0j)

    Sdia.set(0,0, 1.0+0.0j);     Sdia.set(0,1, ex );
    Sdia.set(1,0, ex);           Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );  d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);             d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        d = -(x-0.5*x0)*ex
        dc1_dia[i].set(0,0, 0.0+0.0j);  dc1_dia[i].set(0,1, d);
        dc1_dia[i].set(1,0, d);         dc1_dia[i].set(1,1, 0.0+0.0j);

    

def model4(q, params):
    """

    2-level system with periodic anharmonic potentials 
    

              k*cos(w*x)         V
    Hdia =       V        k*sin(w*x) + D

    Sdia =  I

    Ddia  = 0.0

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha]
            * **params["w"]** ( double ): frequency  [ default: 0.1, units: Bohr^-1]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [ ] 
    default_params = {"k":0.01, "D":0.0, "V":0.005, "w":0.1 }
    comn.check_input(params, default_params, critical_params)
    k,D,V,w = params["k"], params["D"], params["V"], params["w"]

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  
    x = q.get(0)

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*math.cos(x*w)*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));                  Hdia.set(1,1, k*math.sin(x*w)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0,-w*k*math.sin(x*w)*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);                        d1ham_dia[i].set(1,1, w*k*math.cos(x*w)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model4a(Hdia, Sdia, d1ham_dia, dc1_dia, q, params):
    """

    Same as ::funct:```model4``` just different interface
    
              k*cos(w*x)         V
    Hdia =       V        k*sin(w*x) + D

    Sdia =  I

    Ddia  = 0.0

    Args: 
        Hdia ( CMATRIX(2,2) ): diabatic Hamiltonian - updated by this function
        Sdia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states - updated by this function [ identity ] 
        d1ham_dia ( list of 1 CMATRIX(2,2) objects ): derivatives of the diabatic Hamiltonian w.r.t. 
            the nuclear coordinate - updated by this function
        dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis - updated 
            by this function [ zero ]
        q ( double ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["k"]** ( double ): force constante [ default: 0.01, units: Ha]
            * **params["w"]** ( double ): frequency  [ default: 0.1, units: Bohr^-1]
            * **params["D"]** ( double ): gap between the minima of the states 1 and 0, negative 
                value means the state 1 is lower in energy than state 0  [ default: 0.0, units: Ha]
            * **params["V"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]

    Returns:       
        None

    """
    critical_params = [ ] 
    default_params = {"k":0.01, "D":0.0, "V":0.005, "w":0.1 }
    comn.check_input(params, default_params, critical_params)
    k,D,V,w = params["k"], params["D"], params["V"], params["w"]

    x = q

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*math.cos(x*w)*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));                  Hdia.set(1,1, k*math.sin(x*w)*(1.0+0.0j) + D);


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0,-w*k*math.sin(x*w)*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);                        d1ham_dia[i].set(1,1, w*k*math.cos(x*w)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);



def model5(q, params):
    """

    Symmetric Double Well Potential in 2D
    Hdia = A*[0.25*(q1^4 + q2^4) - 0.5*(q1^2 + q2^2) ]
    Sdia = 1.0
    Ddia = 0.0

    Args: 
        q ( MATRIX(2,1) ): coordinates of the particle, ndof = 2
        params ( dictionary ): model parameters

            * **params["A"]** ( double ): scaling parameter to determine the depth
            of the potential well/barrier  [ default: 1.0, units: None]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(1,1) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(1,1) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 2 CMATRIX(1,1) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(1,1) objects ): derivative coupling in the diabatic basis [ zero ] 

    """

    # Hdia and Sdia are ndia x ndia in dimension
    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)

    # d1ham and dc1_dia are ndia x ndia in dimension, but we have nnucl of them
    d1ham_dia = CMATRIXList();
    dc1_dia = CMATRIXList();

    for i in range(0,2):
        d1ham_dia.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) )

    x = q.get(0)
    y = q.get(1)
    x2 = x*x
    y2 = y*y

    Hdia.set(0,0, (0.25*(x2*x2 + y2*y2) - 0.5*(x2 + y2))*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, x*(x2 - 1.0)*(1.0+0.0j) )
    d1ham_dia[1].set(0,0, y*(y2 - 1.0)*(1.0+0.0j) )

    #  <dia| d/dR_0| dia >
    dc1_dia[0].set(0,0, 0.0+0.0j)
    dc1_dia[1].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def model6(q, params):
    """

    2-level system with periodic anharmonic potentials. More general than model 4
    
              A0*cos(w0*x+delta0) + B0         V
    Hdia =               V        A1*cos(w1*x+delta1) + B1

    Sdia =  I

    Ddia  = 0.0

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A0"]** ( double ): amplitude of the state 0 [ default: 0.01, units: Ha]
            * **params["A1"]** ( double ): amplitude of the state 1 [ default: 0.01, units: Ha]
            * **params["w0"]** ( double ): frequency of the potential 0 [ default: 0.1, units: Bohr^-1]
            * **params["w1"]** ( double ): frequency of the potential 1 [ default: 0.1, units: Bohr^-1]
            * **params["delta0"]** ( double ): phase shift of potential 0 [ default: 0.0, units: None]
            * **params["delta1"]** ( double ): phase shift of potential 1 [ default: 0.0, units: None]
            * **params["B0"]** ( double ): energy shift of potential 0 [ default: 0.0, units: Ha]
            * **params["B1"]** ( double ): energy shift of potential 1 [ default: 0.0, units: Ha]
            * **params["V01"]** ( double ): electronic coupling between these diabats [ default: 0.005, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [ ] 
    default_params = {"A0":0.01, "w0":0.1, "delta0":0.0, "B0":0.0,
                      "A1":0.01, "w1":0.1, "delta1":0.0, "B1":0.0,
                      "V01":0.005  }
    comn.check_input(params, default_params, critical_params)

    A0, A1 = params["A0"], params["A1"]
    B0, B1 = params["B0"], params["B1"]
    w0, w1 = params["w0"], params["w1"]
    delta0, delta1 = params["delta0"], params["delta1"]
    V01 = params["V01"]



    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, (A0*math.cos(w0*x+delta0) + B0)*(1.0+0.0j) );   Hdia.set(0,1, V01*(1.0+0.0j));
    Hdia.set(1,0, V01*(1.0+0.0j));                  Hdia.set(1,1, (A1*math.cos(w1*x+delta1) + B1)*(1.0+0.0j));


    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0,-w0*A0*math.sin(w0*x+delta0)*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);                        d1ham_dia[i].set(1,1, -w1*A1*math.sin(w1*x+delta1)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj




def model7(q, params):
    """

    3-level system with periodic anharmonic potentials. 
    
              A0*cos(w0*x+delta0) + B0         V01              V02
    Hdia =               V01        A1*cos(w1*x+delta1) + B1    V12
                         V02                   V12          A2*cos(w2*x+delta2) + B2

    Sdia =  I

    Ddia  = 0.0

    Args: 
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["A0"]** ( double ): amplitude of the state 0 [ default: 0.01, units: Ha]
            * **params["A1"]** ( double ): amplitude of the state 1 [ default: 0.01, units: Ha]
            * **params["A2"]** ( double ): amplitude of the state 2 [ default: 0.01, units: Ha]
            * **params["w0"]** ( double ): frequency of the potential 0 [ default: 0.1, units: Bohr^-1]
            * **params["w1"]** ( double ): frequency of the potential 1 [ default: 0.1, units: Bohr^-1]
            * **params["w2"]** ( double ): frequency of the potential 2 [ default: 0.1, units: Bohr^-1]
            * **params["delta0"]** ( double ): phase shift of potential 0 [ default: 0.0, units: None]
            * **params["delta1"]** ( double ): phase shift of potential 1 [ default: 0.0, units: None]
            * **params["delta2"]** ( double ): phase shift of potential 2 [ default: 0.0, units: None]
            * **params["B0"]** ( double ): energy shift of potential 0 [ default: 0.0, units: Ha]
            * **params["B1"]** ( double ): energy shift of potential 1 [ default: 0.0, units: Ha]
            * **params["B2"]** ( double ): energy shift of potential 2 [ default: 0.0, units: Ha]
            * **params["V01"]** ( double ): electronic coupling between diabats 0 and 1 [ default: 0.005, units: Ha]
            * **params["V02"]** ( double ): electronic coupling between diabats 0 and 2 [ default: 0.005, units: Ha]
            * **params["V12"]** ( double ): electronic coupling between diabats 1 and 2 [ default: 0.005, units: Ha]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(3,3) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(3,3) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(3,3) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(3,3) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = [ ] 
    default_params = {"A0":0.01, "w0":0.1, "delta0":0.0, "B0":0.0,
                      "A1":0.01, "w1":0.1, "delta1":0.0, "B1":0.0,
                      "A2":0.01, "w2":0.1, "delta2":0.0, "B2":0.0,
                      "V01":0.005,  "V02":0.005, "V12":0.005  }
    comn.check_input(params, default_params, critical_params)

    A0, A1, A2 = params["A0"], params["A1"], params["A2"]
    B0, B1, B2 = params["B0"], params["B1"], params["B2"]
    w0, w1, w2 = params["w0"], params["w1"], params["w2"]
    delta0, delta1, delta2 = params["delta0"], params["delta1"], params["delta2"]
    V01, V02, V12 = params["V01"], params["V02"], params["V12"]


    Hdia = CMATRIX(3,3)
    Sdia = CMATRIX(3,3)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(3,3))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(3,3))
  

    x = q.get(0)

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);  Sdia.set(0,2, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);  Sdia.set(1,2, 0.0+0.0j);
    Sdia.set(2,0, 0.0+0.0j);  Sdia.set(2,1, 0.0+0.0j);  Sdia.set(2,2, 1.0+0.0j);

    Hdia.set(0,0, (A0*math.cos(w0*x+delta0) + B0)*(1.0+0.0j) ); 
    Hdia.set(0,1, V01*(1.0+0.0j));
    Hdia.set(0,2, V02*(1.0+0.0j));

    Hdia.set(1,0, V01*(1.0+0.0j));
    Hdia.set(1,1, (A1*math.cos(w1*x+delta1) + B1)*(1.0+0.0j));
    Hdia.set(1,2, V12*(1.0+0.0j));

    Hdia.set(2,0, V02*(1.0+0.0j));
    Hdia.set(2,1, V12*(1.0+0.0j));
    Hdia.set(2,2, (A2*math.cos(w2*x+delta2) + B2)*(1.0+0.0j));



    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0,-w0*A0*math.sin(w0*x+delta0)*(1.0+0.0j) );
        d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(0,2, 0.0+0.0j);

        d1ham_dia[i].set(1,0, 0.0+0.0j); 
        d1ham_dia[i].set(1,1, -w1*A1*math.sin(w1*x+delta1)*(1.0+0.0j));
        d1ham_dia[i].set(1,2, 0.0+0.0j); 

        d1ham_dia[i].set(2,0, 0.0+0.0j); 
        d1ham_dia[i].set(2,1, 0.0+0.0j); 
        d1ham_dia[i].set(2,2, -w2*A2*math.sin(w2*x+delta2)*(1.0+0.0j));


        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);  dc1_dia[i].set(0,2, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);  dc1_dia[i].set(1,2, 0.0+0.0j);
        dc1_dia[i].set(2,0, 0.0+0.0j);   dc1_dia[i].set(2,1, 0.0+0.0j);  dc1_dia[i].set(2,2, 0.0+0.0j);


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


