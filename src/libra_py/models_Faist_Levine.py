#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file models_Faist_Levine.py 
#
# This module implements the 1D, 2-level models of Faist and Levine
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



def Faist_Levine(q, params):
    """
    Faist and Levine, 2-level, 1-dim. potential describing the
    collision of alkali and halogen atoms

    \param[in] q [1 x 1, MATRIX] coordinates of the particles 
    \param[in] params [dictionary] parameters of the model. 

    Should contain the key - value pairs: 
        key          value        description

      "A_cov"    -    (double)      
      "A_ion"    -    (double)
      "B_cov"    -    (double)      
      "B_ion"    -    (double)
      "C_cov"    -    (double)      
      "C_ion"    -    (double)
      "rho_cov"  -    (double)      
      "rho_ion"  -    (double)
      "alp_M+"   -    (double)
      "alp_X-"   -    (double)
      "E_th"     -    (double)
      "A"        -    (double)
      "rho"      -    (double)

    Ref: Faist, M. B.; Levine, R. D. JCP 1976, 64, 2953
    
    """

    ndof = q.num_of_rows  # the number of nuclear DOFs

    if ndof != 1:
        print "Error: The Faist-Levine Hamiltonian takes 1 nuclear DOFs only. Given = ", ndof
        print "Exiting..."
        sys.exit(0)

    A_cov = params["A_cov"]
    A_ion = params["A_ion"]
    B_cov = params["B_cov"]
    B_ion = params["B_ion"]
    C_cov = params["C_cov"]
    C_ion = params["C_ion"]
    rho_cov = params["rho_cov"]
    rho_ion = params["rho_ion"]
    alp_Mp = params["alp_M+"]
    alp_Xm = params["alp_X-"]
    E_th = params["E_th"]
    A = params["A"]
    rho = params["rho"]


    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in xrange(1):
        obj.d1ham_dia.append( CMATRIX(2,2) )
        obj.dc1_dia.append( CMATRIX(2,2) )


    #=========== Energies & Derivatives ===============
    R = q.get(0)

    p1 = 1.0/R
    p2 = p1 * p1
    p4 = p2 * p2
    p6 = p4 * p2
    p8 = p6 * p2
    p12 = p6 * p6


    # H_00
    e = math.exp(-R/rho_cov)
    b12 = math.pow(B_cov)
    pb12 = b12 * p12

    z = (A_cov + pb12 ) * e - C_cov * p8
    dzdx =  (-12.0 * pb12 / R ) * e
    dzdx += (-1.0/rho_cov) * (A_cov + pb12 ) * e
    dzdx += 8.0 * C_cov * p8 / R

    obj.ham_dia.set(0,0, z*(1.0+0.0j))
    obj.d1ham_dia[0].set(0,0, dzdx*(1.0+0.0j))


    # H_11
    e = math.exp(-R/rho_ion)
    b8 = math.pow(B_ion)
    pb8 = b8 * p8

    z = (A_ion + pb8 ) * e - C_ion * p6 - p1 - 0.5*(alp_Mp + alp_Xm)*p4 - 2.0*alp_Mp*alp_Xm*p6/R + E_th
    dzdx =  (-8.0 * pb8 / R ) * e
    dzdx += (-1.0/rho_ion) * (A_ion + pb8 ) * e
    dzdx += 6.0 * C_ion * p6 / R
    dzdx += p2
    dzdx += 2.0*(alp_Mp + alp_Xm)*p4 / R 
    dzdx += 14.0*alp_Mp*alp_Xm*p8

    obj.ham_dia.set(1,1, z * (1.0+0.0j))
    obj.d1ham_dia[0].set(1,1, dzdx * (1.0+0.0j))



    # H_01 and H_10
    z = A*math.exp(-R/rho)
    dzdx = (-1.0/rho) * z

    obj.ham_dia.set(0,1, z * (1.0+0.0j) )
    obj.d1ham_dia[0].set(0,1, dzdx * (1.0+0.0j) )

    obj.ham_dia.set(1,0, z * (1.0+0.0j) )
    obj.d1ham_dia[0].set(1,0, dzdx * (1.0+0.0j) )


    return obj



def get_Faist_Levine_LiI():
    """
    Parameters for Li + I collision
    
    Ref: Faist, M. B.; Levine, R. D. JCP 1976, 64, 2953

    """

    A = 1.889725989  # 1 Angstrom in Bohr
    eV = 0.036749309 # 1 eV in Ha

    params = {}

    params["A_cov"] = 1100.0 * eV
    params["A_ion"] = 1052.0 * eV
    params["B_cov"] = 2.996 * math.pow(eV, 1.0/12.0) * A
    params["B_ion"] = 1.839 * math.pow(eV, 1.0/8.0) * A
    params["C_cov"] = 1191.2 * eV * math.pow(A, 8)
    params["C_ion"] = 0.823 * eV * math.pow(A, 6)
    params["rho_cov"] = 0.44 * eV
    params["rho_ion"] = 0.3786 * eV
    params["alp_M+"] = 0.029 * math.pow(A, 3)
    params["alp_X-"] = 6.431 * math.pow(A, 3)
    params["E_th"] = 2.326 * eV
    params["A"] = 17.08 * eV
    params["rho"] = 1.239 * A  ## I'm not sure why they give it in A^-1


    return params
   


def get_Faist_Levine_NaI():
    """
    Parameters for Na + I collision
    
    Ref: Faist, M. B.; Levine, R. D. JCP 1976, 64, 2953

    """

    A = 1.889725989  # 1 Angstrom in Bohr
    eV = 0.036749309 # 1 eV in Ha

    params = {}

    params["A_cov"] = 3150.0 * eV
    params["A_ion"] = 2760.0 * eV
    params["B_cov"] = 2.647 * math.pow(eV, 1.0/12.0) * A
    params["B_ion"] = 2.398 * math.pow(eV, 1.0/8.0) * A
    params["C_cov"] = 1000.0 * eV * math.pow(A, 8)
    params["C_ion"] = 11.3 * eV * math.pow(A, 6)
    params["rho_cov"] = 0.435 * eV
    params["rho_ion"] = 0.3489 * eV
    params["alp_M+"] = 0.408 * math.pow(A, 3)
    params["alp_X-"] = 6.431 * math.pow(A, 3)
    params["E_th"] = 2.075 * eV
    params["A"] = 17.08 * eV
    params["rho"] = 1.239 * A  ## I'm not sure why they give it in A^-1


    return params
   
