# *********************************************************************************
# * Copyright (C) 2020 Story Temen, Alexey V. Akimov
# * Copyright (C) 2018-2019 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 2 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: models_Holstein
   :platform: Unix, Windows
   :synopsis: This module implements the Henon-Heiles Hamiltonians
.. moduleauthor:: Alexey V. Akimov, Story Temen

"""

import os
import sys
import math
import copy
import torch

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units


class tmp:
    pass

def Holstein2(q, params):
    """
    n-state model

    H_nn = E_n + 0.5*k*(x-x_n)^2

    H_n,n+1 = H_n+1,n = V,
    H_n,m = 0, otherwise

    Args:
        q ( Tensor(nbeads, nnucl) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["E_n"]** ( list of doubles ):  [ default: [0.0, 0.001, 0.001, 0.001], units: Ha]
            * **params["x_n"]** ( list of doubles ):  [ default: [0.0, 1.0, 1.0, 1.0], units: Bohr]
            * **params["k_n"]** ( list of doubles ):  [ default: [0.001, 0.001, 0.001, 0.001], units: Ha/Bohr^2]
            * **params["V"]**   ( double ):  [ default: 0.001, units: Ha]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(4,4) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(4,4) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(4,4) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(4,4) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = []
    default_params = {"V": 0.001, "E_n":[0.0, 0.001, 0.001], "x_n":[0.0, 1.0, 2.0], "k_n":[0.001, 0.001, 0.001] }
    comn.check_input(params, default_params, critical_params)

    E_n = params["E_n"]
    x_n = params["x_n"]
    k_n = params["k_n"]
    V = params["V"]

    nbeads, nnucl = q.shape[0], q.shape[1]
    nstates = len(E_n)

    #x = q.col(indx).get(0)

    obj = tmp()
    obj.ham_dia = torch.zeros(nbeads,nstates,nstates,dtype=torch.complex128)
    obj.ovlp_dia = torch.eye(nstates, dtype=torch.complex128).unsqueeze(0).repeat(nbeads, 1, 1) # nbeads of nstates x nstates identity matrices
    obj.d1ham_dia = torch.zeros(nbeads,nnucl, nstates,nstates,dtype=torch.complex128)
    obj.dc1_dia = torch.zeros(nbeads,nnucl, nstates,nstates, dtype=torch.complex128)

    for b in range(nbeads):
        for i in range(nstates):
            obj.ham_dia[b, i, i] = E_n[i] * (1.0 + 0.0j) 

            for n in range(nnucl):
                obj.ham_dia[b, i, i] +=  0.5 * k_n[i] * (q[b, n] - x_n[i])**2 * (1.0 + 0.0j)
                obj.d1ham_dia[b, n, i, i] = k_n[i] * (q[b, n] - x_n[i]) * (1.0 + 0.0j)
        
        for i in range(nstates-1):    
            obj.ham_dia[b, i, i+1] = V * (1.0 + 0.0j)
            obj.ham_dia[b, i+1, i] = V * (1.0 + 0.0j)



    return obj


#x = torch.Tensor(2,2,2).zero_()
#print(x, x.shape)

