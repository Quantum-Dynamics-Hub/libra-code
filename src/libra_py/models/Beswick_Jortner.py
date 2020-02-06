#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file models_Beswick_Jortner.py 
#
# This module implements the Beswick-Jortner potential describing dissociation of ICN molecule
# or collinear interaction of I and CN
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


def Beswick_Jortner(q, params):
    """
    Implementation of the Henon-Heiles potential

    Args:
        q ( MATRIX(ndof, 1) ): coordinate of the particle, ndof = 2 
               These DOFs are as follows:
               q[0] = r - C-N distance
               q[1] = R - the distance between I and the center of mass of CN diatomics
        params ( dictionary ): parameters of the model

    References:
        (1) A. Beswick and J. Jortner, Chern. Phys. 1977, 24, 1 
        (2) Brown, R. C.; Heller, E. J. Classical Trajectory Approach to Photodissociation: The Wigner Method. 
        J. Chem. Phys. 1981, 75, 186–188.
    """

    eV = 0.036749309 #  Ha
    Angst = 1.889725989  # Bohr

    K = 74.4434 * eV / Angst   # force constant
    r0 = 1.2327 * eV           # equilibrium C-N distance
    T0 = 4.36994 * eV  
    A = 200000 * eV
    a = 6.68 * (1/Angst)

    """
    Here are the parameters to determine the effective masses of the DOFs
              Ground state:              Excited state 
    w1 (cm-1)      499                    w(CN, cm^-1)   1763
    w2 (cm-1)     2237
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

    r = q.get(0)
    R = q.get(1)

    Hdia.set(0,0,( T0 + 0.5*k*(r-r0)**2 + A*exp(a* (r*m_c/(m_c + m_n) - R) )  )*(1.0+0.0j) )
    Sdia.set(0,0,1.0+0.0j)

    #  d Hdia / dR_0
    d1ham_dia[0].set(0,0, ( x + lam*(y2 - x2) + 0.25 * lam2*(x2 + y2)*x )*(1.0+0.0j))
    d1ham_dia[1].set(0,0, ( y + lam*2.0*y*x   + 0.25 * lam2*(x2 + y2)*y )*(1.0+0.0j))

    #  <dia| d/dR_0| dia >
    dc1_dia[0].set(0,0, 0.0+0.0j)
    dc1_dia[1].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


