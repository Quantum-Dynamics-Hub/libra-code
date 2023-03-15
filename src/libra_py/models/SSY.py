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
.. module:: models_SSY
   :platform: Unix, Windows
   :synopsis: This module implements the 2D, 2-level model of Shenvi-Subotnik-Yang
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


def SSY(q, params, full_id):
    """

    The Hamiltonian of Shenvi-Subotnik-Yang, 2-level, 2-dim. problem

    Reference: Shenvi, N.; Subotnik, J.; Yang, W. JCP 2011, 135, 024101

    Args:
        q ( MATRIX(2,1) ): coordinates of the particle, ndof = 2 
        params ( dictionary ): model parameters

            * **params["E0"]** ( double ):  [ default: 0.05, units: Ha]
            * **params["A"]** ( double ): [ default: 0.15, units: Ha]
            * **params["B"]** ( double ): [ default: 0.14, units: Bohr^-2]
            * **params["C"]** ( double ): [ default: 0.015, units: Ha]
            * **params["D"]** ( double ): [ default: 0.06, units: Bohr^-2]


    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian 
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 2 CMATRIX(2,2) objects ): 
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = [ ] 
    default_params = {"E0":0.05, "A":0.15, "B":0.14, "C":0.015, "D":0.06 }
    comn.check_input(params, default_params, critical_params)

    E0 = params["E0"]
    A = params["A"]
    B = params["B"]
    C = params["C"]
    D = params["D"]


    ndof = q.num_of_rows  # the number of nuclear DOFs
    if ndof != 2:
        print("Error: The SSY Hamiltonian takes 2 nuclear DOFs only. Given = ", ndof)
        print("Exiting...")
        sys.exit(0)


    obj = tmp()
    obj.ham_dia = CMATRIX(2,2)
    obj.ovlp_dia = CMATRIX(2,2);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in range(0,2):
        obj.d1ham_dia.append( CMATRIX(2,2) )
        obj.dc1_dia.append( CMATRIX(2,2) )


    #=========== Energies & Derivatives ===============

    Id = Cpp2Py(full_id)
    indx = Id[-1]
    x, y = q.get(0, indx), q.get(1, indx)

    # H_00
    obj.ham_dia.set(0,0, E0*(-1.0+0.0j))
    obj.d1ham_dia[0].set(0,0, 0.0+0.0j)
    obj.d1ham_dia[1].set(0,0, 0.0+0.0j)

    # H_11
    z = -A*math.exp(-B * (0.75*(x+y)**2 + 0.25*(x-y)**2) )
    dzdx = -B * (1.5*(x+y) + 0.5*(x-y)) * z
    dzdy = -B * (1.5*(x+y) - 0.5*(x-y)) * z
    obj.ham_dia.set(1,1, z * (1.0+0.0j))
    obj.d1ham_dia[0].set(1,1, dzdx * (1.0+0.0j))
    obj.d1ham_dia[1].set(1,1, dzdy * (1.0+0.0j))


    # H_01 and H_10
    z = C*math.exp(-D * (0.25*(x+y)**2 + 0.75*(x-y)**2) )
    dzdx = -D * (0.5*(x+y) + 1.5*(x-y)) * z
    dzdy = -D * (0.5*(x+y) - 1.5*(x-y)) * z

    obj.ham_dia.set(0,1, z * (1.0+0.0j) )
    obj.d1ham_dia[0].set(0,1, dzdx * (1.0+0.0j) )
    obj.d1ham_dia[1].set(0,1, dzdy * (1.0+0.0j) )

    obj.ham_dia.set(1,0, z * (1.0+0.0j) )
    obj.d1ham_dia[0].set(1,0, dzdx * (1.0+0.0j) )
    obj.d1ham_dia[1].set(1,0, dzdy * (1.0+0.0j) )


    return obj

