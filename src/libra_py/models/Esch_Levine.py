#*********************************************************************************
#* Copyright (C) 2021 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#***********************************************************************************
"""
.. module:: models_Esch_Levine
   :platform: Unix, Windows
   :synopsis: This module implements the Hamiltonians of Esch and LEvine
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



def get_JCP_2020_params(set_indx):
    """
    Parameters from:
    Esch, M. P.; Levine, B. G. J. Chem. Phys. 153, 114104 (2020); doi: 10.1063/5.0022529

    Args:
        None

    Returns:
        dictionary: params, will contain the parameters:
            * **params["w0"]** ( double ): slope of the negative diabat [ default: 0.25, units: Ha]
            * **params["w1"]** ( double ): slope of the positive diabats [ default: 0.025, units: Ha]
            * **params["V"]**  ( double ): diabatic coupling [ default: 0.005, units: Ha]
            * **params["eps"]** ( double ): the amount of the energy shift for the group of states [ default: 0.0, units: Ha]
            * **params["i_crit"]** ( int ): the index of the level starting from which the group of states is shifted [ default: 2 ]
            * **params["nstates"]** ( int ): the number of states, M [ default : 2]

    """
    params = {}

    if set_indx==1:
        params = { "w0":0.25, "w1":0.025, "V":0.005, "eps":0.0, "i_crit":0, "nstates":9, "delta":0.01 }
    elif set_indx==2:
        params = { "w0":0.25, "w1":0.025, "V":0.005, "eps":0.08, "i_crit":5, "nstates":9, "delta":0.01 }
    elif set_indx==3:
        params = { "w0":0.25, "w1":0.025, "V":0.005, "eps":0.0, "i_crit":0, "nstates":17, "delta":0.005 }
    elif set_indx==4:
        params = { "w0":0.25, "w1":0.025, "V":0.005, "eps":0.0, "i_crit":0, "nstates":17, "delta":0.0025 }
    elif set_indx==5:
        params = { "w0":0.25, "w1":0.025, "V":0.005, "eps":0.0, "i_crit":0, "nstates":17, "delta":0.00125 }
    elif set_indx==6:
        params = { "w0":0.25, "w1":0.025, "V":0.005, "eps":0.0, "i_crit":0, "nstates":17, "delta":0.000625 }



    return params



def JCP_2020(q, params, full_id):
    """
    The models from 
    Esch, M. P.; Levine, B. G. J. Chem. Phys. 153, 114104 (2020); doi: 10.1063/5.0022529

    M-state model

    H_00 = -w0 * x,  w0 >0
    H_ii = w1 * x  - i * delta - eps_i ,  0 < i <= M-1,  w1 > 0, delta > 0
    where eps_i = eps, for i >= i_crit
    
    H_0,n = H_n,0 = V,  0 < n <= M-1
    H_n,m = 0, otherwise

    Args:
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["w0"]** ( double ): slope of the negative diabat [ default: 0.25, units: Ha]
            * **params["w1"]** ( double ): slope of the positive diabats [ default: 0.025, units: Ha]
            * **params["delta"]** ( double ): energy shift between nearby diabats [ default: 0.01, units: Ha]
            * **params["V"]**  ( double ): diabatic coupling [ default: 0.005, units: Ha]
            * **params["eps"]** ( double ): the amount of the energy shift for the group of states [ default: 0.0, units: Ha]
            * **params["i_crit"]** ( int ): the index of the level starting from which the group of states is shifted [ default: 2 ]
            * **params["nstates"]** ( int ): the number of states, M [ default : 2]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(M,M) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(M,M) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(M,M) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(M,M) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = [ ]
    default_params = { "w0":0.25, "w1":0.025, "V":0.005, "eps":0.0, "i_crit":2, "nstates":2, "delta":0.01 }
    comn.check_input(params, default_params, critical_params)

    w0 = params["w0"]
    w1 = params["w1"]
    delta = params["delta"]
    V = params["V"]
    n = params["nstates"]
    eps = params["eps"]
    i_crit = params["i_crit"]

    Hdia = CMATRIX(n,n)
    Sdia = CMATRIX(n,n)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(n,n) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(n,n) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    x = q.col(indx).get(0)

    Sdia.identity()

    Hdia.set(0,0,  -w0 * x * (1.0+0.0j) )
    for i in range(1,n):
        shift = 0.0
        if i>=i_crit:
            shift = eps
        Hdia.set(i,i,  (w1 * x - i * delta - shift) * (1.0+0.0j) )

    for i in range(1,n):
        Hdia.set(0,i,  V * (1.0+0.0j) )
        Hdia.set(i,0,  V * (1.0+0.0j) )

    for k in [0]:
        #  d Hdia / dR_0
        d1ham_dia[k].set(0,0, -w0 * (1.0+0.0j) )
        for i in range(1, n):
            d1ham_dia[k].set(i,i, w1*(1.0+0.0j) )


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
