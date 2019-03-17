#*********************************************************************************                     
#* Copyright (C) 2018 Brendan A. Smith, Alexey V. Akimov                                                   
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
    The first of the two model potential described by L. Wang, C.C. Martens, 
    and Y. Zheng, "Entangled trajectory molecular dynamics in 
    multidimensional systems: Two-dimensional quantum tunneling through 
    the Eckart barrier" J. Chem. Phys. 137, 34113 (2012)

    This potential has seperable dof

    # Define system specific parameters 
    Va = barrier height     = 0.00625    
    Vb = harmonic potential = 0.0106

    Sdia = 1.0
    Ddia = 0.0
    """

    # Define potential specific constants
    Va = 0.00625
    Vb = 0.0106

    # Hdia and Sdia are ndia x ndia in dimension
    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)

    # d1ham and dc1_dia are ndia x ndia in dimension, but we have nnucl of them
    d1ham_dia = CMATRIXList();
    dc1_dia = CMATRIXList();

    for i in xrange(2):
        d1ham_dia.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) )

    x = q.get(0)
    y = q.get(1)
    x2 = x*x
    y2 = y*y

    # z = sech(2x)
    # z2 = sech^2(2x)    
    z = (2.0 * math.cosh(2.0*x))/(1.0 + math.cosh(4.0*x))
    z2 = z*z

    Hdia.set(0,0,( Va*z2 + 0.5*Vb*y*y )*(1.0+0.0j))
    Sdia.set(0,0,1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, ( ( -4.0*Va*math.tanh(2.0*x)*z2 ) )*(1.0+0.0j))
    d1ham_dia[1].set(0,0, ( Vb*y )*(1.0+0.0j))

    #  <dia| d/dR_0| dia >
    dc1_dia[0].set(0,0, 0.0+0.0j)
    dc1_dia[1].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj

def model2(q, params):
    """
    The second of the two model potential described by L. Wang, C.C. Martens, 
    and Y. Zheng, "Entangled trajectory molecular dynamics in 
    multidimensional systems: Two-dimensional quantum tunneling through 
    the Eckart barrier" J. Chem. Phys. 137, 34113 (2012)

    Specifically, this potential is Model II 

    Hdia = Va*sech^2*(2q1) + 0.5*Vb*[q2 - Vc(q1^2 - 1)]^2

    Va = barrier height     = 0.00625    
    Vb = harmonic potential = 0.0106
    Vc = coupling constant  = 0.4  

    Sdia = 1.0
    Ddia = 0.0
    """

    # Define potential specific constants
    Va = 0.00625    
    Vb = 0.0106
    Vc = 0.4    

    # Hdia and Sdia are ndia x ndia in dimension
    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)

    # d1ham and dc1_dia are ndia x ndia in dimension, but we have nnucl of them
    d1ham_dia = CMATRIXList();
    dc1_dia = CMATRIXList();

    for i in xrange(2):
        d1ham_dia.append( CMATRIX(1,1) )
        dc1_dia.append( CMATRIX(1,1) )

    x = q.get(0)
    y = q.get(1)
    x2 = x*x
    y2 = y*y

    # z = sech(2x)
    # z2 = sech^2(2x)    
    z = (2.0 * math.cosh(2.0*x))/(1.0 + math.cosh(4.0*x))
    z2 = z*z

    Hdia.set(0,0,( Va*z2 + 0.5*Vb*(y+Vc*(x*x-1.0))**2 )*(1.0+0.0j))
    Sdia.set(0,0,1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, ( ( -4.0*Va*math.tanh(2.0*x)*z2 ) + ( 2.0*Vb*Vc*x*(y + Vc*x*x - Vc) ) )*(1.0+0.0j))
    d1ham_dia[1].set(0,0, ( Vb * ( y + Vc*(x*x-1.0) ) ) * (1.0+0.0j))

    #  <dia| d/dR_0| dia >
    dc1_dia[0].set(0,0, 0.0+0.0j)
    dc1_dia[1].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
