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
.. module:: Morse
   :platform: Unix, Windows
   :synopsis: This module implements various model Hamiltonians composed of N displaced
       and scaled Morse potentials
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

def set_Coronado_Xing_Miller_params(params, model_indx):
    """

    This function will re-set the corresponding parameters of the 
    `params` dictionary according to the model Hamiltonian select by
    `model_indx`

    The 3-state Morse potential described in:

    Corondao, E. A.; Xing, J.; Miller, W. H. Chem. Phys. Lett. 2001, 349, 521-529

    Note: the definition of `alpha` and `beta` variables chosen in this module
    is opposite to the one used in the paper

    """

    if model_indx == 1:
        tmp  = { "D": [0.003, 0.004, 0.003],  "alpha": [0.65, 0.6, 0.65],   "x_n": [5.00, 4.00, 6.00],    "E": [0.00, 0.01, 0.006],
                 "V": [ [0.000, 0.002, 0.000], 
                        [0.002, 0.000, 0.002],
                        [0.000, 0.002, 0.000] ],
                 "x_nm": [ [0.00, 3.40, 0.00], 
                           [3.40, 0.00, 4.80],
                           [0.00, 4.80, 0.00] ],
                 "beta": [ [0.00, 16.00, 0.00], 
                           [16.00, 0.00, 16.00],
                           [0.00, 16.00, 0.00] ]
                }

    elif model_indx == 2:
        tmp  = { "D": [0.02, 0.01, 0.003],  "alpha": [0.65, 0.4, 0.65],   "x_n": [4.50, 4.00, 4.40],    "E": [0.00, 0.01, 0.02],
                 "V": [ [0.000, 0.005, 0.005], 
                        [0.005, 0.000, 0.000],
                        [0.005, 0.000, 0.000] ],
                 "x_nm": [ [0.00, 3.66, 3.34], 
                           [3.66, 0.00, 0.00],
                           [3.34, 0.00, 0.00] ],
                 "beta": [ [0.00, 32.00, 32.00], 
                           [32.00, 0.00, 0.00],
                           [32.00, 0.00, 0.00] ]
                }
          
    elif model_indx == 3:
        tmp  = { "D": [0.02, 0.02, 0.003],  "alpha": [0.4, 0.65, 0.65],   "x_n": [4.0, 4.50, 6.0],    "E": [0.02, 0.00, 0.02],
                 "V": [ [0.000, 0.005, 0.005], 
                        [0.005, 0.000, 0.000],
                        [0.005, 0.000, 0.000] ],
                 "x_nm": [ [0.00, 3.4, 4.97], 
                           [3.4, 0.00, 0.00],
                           [4.97, 0.00, 0.00] ],
                 "beta": [ [0.00, 32.00, 32.00], 
                           [32.00, 0.00, 0.00],
                           [32.00, 0.00, 0.00] ]
                }


    params.update(tmp)    

 

def general(q, params, full_id):
    """
    H_nn = D * [1 - exp(-alpha_n * (x-x_n))]^2 + E
    H_nm = V * exp(-beta_nm * (x-x_nm)^2 )

    Args:
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["D"]** ( list of doubles ):  [ default: [0.001, 0.001, 0.001], units: Ha]
            * **params["alpha"]** ( list of double ):  [ default:  [0.0, 0.0, 0.0], units: 1/Bohr^2]
            * **params["E"]** ( list of doubles ):  [ default: [0.0, -0.001, -0.002], units: Ha]
            * **params["x_n"]** ( list of doubles ):  [ default: [0.0, 1.0, 2.0], units: Bohr]

            * **params["V"]**   ( list of lists of double ):  [ default:  [[0.001]*3]*3, units: Ha]
            * **params["beta"]** ( list of lists of double ):  [ default:  [[0.0]*3]*3, units: 1/Bohr^2]
            * **params["x_nm"]**  ( list of lists of double ):  [ default:  [[0.0]*3]*3, units: Bohr]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(3,3) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(3,3) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(3,3) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.d2ham_dia ( list of 1 CMATRIX(3,3) objects ):
                second derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(3,3) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = ["E"]
    default_params = { "D": [0.00, 0.00, 0.00],
                       "alpha": [0.00, 0.00, 0.00],
                       "x_n": [0.00, 0.00, 0.00],

                       "V": [ [0.001, 0.001, 0.001], 
                              [0.001, 0.001, 0.001],
                              [0.001, 0.001, 0.001] ],
                       "beta": [ [0.00, 0.00, 0.00], 
                                  [0.00, 0.00, 0.00],
                                  [0.00, 0.00, 0.00] ],
                       "x_nm": [ [0.00, 0.00, 0.00], 
                                 [0.00, 0.00, 0.00],
                                 [0.00, 0.00, 0.00] ],
                     }
    comn.check_input(params, default_params, critical_params)

    E = params["E"]
    D = params["D"]
    alpha = params["alpha"]
    x_n = params["x_n"]
    V = params["V"]
    beta = params["beta"]
    x_nm = params["x_nm"]

    n = len(E)

    Hdia = CMATRIX(n,n)
    Sdia = CMATRIX(n,n)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(n,n) )
    d2ham_dia = CMATRIXList();  d2ham_dia.append( CMATRIX(n,n) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(n,n) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    x = q.col(indx).get(0)

    Sdia.identity()

    # Diagonals
    for i in range(n):
        bo = math.exp(-alpha[i]*(x-x_n[i]))
        Hdia.set(i,i,  (E[i] + D[i] * (1.0 - bo)**2 ) * (1.0+0.0j) )

        for k in [0]:
            #  d Hdia / dR_0
            d1ham_dia[k].set(i,i, 2.0*alpha[i]*D[i]*bo*(1.0 - bo) * (1.0+0.0j) )


    # Off-diagonals
    for i in range(n):
        for j in range(n):           
            if i!=j:
                d = V[i][j] * math.exp(-beta[i][j] * (x-x_nm[i][j])**2)

                Hdia.set(i,j,  d * (1.0+0.0j) )

                for k in [0]:
                    d1ham_dia[k].set(i,j,  -2.0*beta[i][j] * (x-x_nm[i][j]) * d * (1.0+0.0j) )



    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.d2ham_dia = d2ham_dia
    obj.dc1_dia = dc1_dia

    return obj

