#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file Libra.py 
#
#  Original Libra models - well, they are not necessarily introduced here for the 
#  first time, but these are the models I define for the internal test purposes
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



def model1(q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    """


    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2) )
  

    x = q.get(0)
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

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model1a(Hdia, Sdia, d1ham_dia, dc1_dia, q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

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

    

def model2(q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia = 0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)
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
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1,-0.1+0.0j);
        dc1_dia[i].set(1,0, 0.1+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model2a(Hdia, Sdia, d1ham_dia, dc1_dia, q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia = 0.0

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
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1,-0.1+0.0j);
        dc1_dia[i].set(1,0, 0.1+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);



def model3(q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =    1                           0.05*exp(-(x-0.5*x0)^2)
             0.05*exp(-(x-0.5*x0)^2)               1

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia !=0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)
    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    ex = math.exp(-(x-0.5*x0)**2)*(0.05+0.0j)

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


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model3a(Hdia, Sdia, d1ham_dia, dc1_dia, q, params):
    """
              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =    1                           0.05*exp(-(x-0.5*x0)^2)
             0.05*exp(-(x-0.5*x0)^2)               1

    Ddia != 0.0, but Ddia + Ddia.H() = dSdia/dR, with dSdia !=0.0

    """


    x = q
    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    ex = math.exp(-(x-0.5*x0)**2)*(0.05+0.0j)

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
              k*cos(w*x)         V
    Hdia =       V        k*sin(w*x) + D

    Sdia =  I

    Ddia  = 0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))
  

    x = q.get(0)
    k, w, V = params["k"], params["omega"], params["V"]

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
              k*cos(w*x)         V
    Hdia =       V        k*sin(w*x) + D

    Sdia =  I

    Ddia  = 0.0

    """

    x = q
    k, w, V = params["k"], params["omega"], params["V"]

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


