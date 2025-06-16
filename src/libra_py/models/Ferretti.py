# *********************************************************************************
# * Copyright (C) 2023 Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: Ferretti
   :platform: Unix, Windows
   :synopsis: This module implements Ferretti's 2D model: Ferretti, A.; Granucci, G.; Lami, A.; Persico, M.; Villani, G.
       Quantum Mechanical and Semiclassical Dynamics at a Conical Intersection. The Journal of Chemical Physics 1996, 104,
       5517–5527. https://doi.org/10.1063/1.471791
.. moduleauthor:: Alexey V. Akimov

"""

import os
import sys
import math
import copy

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units


class tmp:
    pass


def Ferretti(q, params, full_id=None):
    """
    2-level, 2D model of Ferretti et al: Ferretti, A.; Granucci, G.; Lami, A.; Persico, M.; Villani, G.
       Quantum Mechanical and Semiclassical Dynamics at a Conical Intersection. The Journal of Chemical Physics 1996, 104,
       5517–5527. https://doi.org/10.1063/1.471791

    H_11(X,Y) = 0.5 * Kx * (X-X1)**2 + 0.5 * Ky * Y**2
    H_22(X,Y) = 0.5 * Kx * (X-X2)**2 + 0.5 * Ky * Y**2 + Delta
    H_12(X,Y) = gamma * Y * exp(-alpha * (X-X3)**2  - beta * Y**2 )
    Kx = 0.02   Ky = 0.10  Delta = 0.01
    X1 = 4   X2 = X3 = 3
    alpha = 3.0     beta = 1.5    gamma = 0.005 to 0.08

    Use it with masses: Mx = 20 000,  My = 6667

    Args:
        q ( MATRIX(2, 1) ): coordinates of the classical particles
        params ( dictionary ): model parameters, should contain:

            * **params["X1"]** ( double ):[ units: Bohr, default: 4 ]
            * **params["X2"]** ( double ):[ units: Bohr, default: 3 ]
            * **params["X3"]** ( double ):[ units: Bohr, default: 3 ]
            * **params["Kx"]** ( double ):[ units: Ha/Bohr^2, default: 0.02 ]
            * **params["Ky"]** ( double ):[ units: Ha/Bohr^2, default: 0.10 ]
            * **params["Delta"]** ( double ):[ units: Ha, default: 0.01 ]
            * **params["alpha"]** ( double ):[ units: Bohr^-2, default: 3 ]
            * **params["beta"]** ( double ):[ units: Bohr^-2, default: 1.5 ]
            * **params["gamma"]** ( double ):[ units: Ha, default: 0.005 ]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 2 CMATRIX(2,2) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]


    """

    # Define potential specific constants
    critical_params = []
    default_params = {
        "X1": 4.0,
        "X2": 3.0,
        "X3": 3.0,
        "Kx": 0.02,
        "Ky": 0.1,
        "Delta": 0.01,
        "alpha": 3.0,
        "beta": 1.5,
        "gamma": 0.005}
    comn.check_input(params, default_params, critical_params)

    X1 = params["X1"]
    X2 = params["X2"]
    X3 = params["X3"]
    Kx = params["Kx"]
    Ky = params["Ky"]
    alpha = params["alpha"]
    beta = params["beta"]
    gamma = params["gamma"]
    Delta = params["Delta"]

    ndof = q.num_of_rows  # the number of nuclear DOFs

    obj = tmp()
    obj.ham_dia = CMATRIX(2, 2)
    obj.ovlp_dia = CMATRIX(2, 2)
    obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()
    obj.d2ham_dia = CMATRIXList()

    for i in range(0, ndof):
        obj.d1ham_dia.append(CMATRIX(2, 2))
        obj.dc1_dia.append(CMATRIX(2, 2))

        for j in range(0, ndof):
            obj.d2ham_dia.append(CMATRIX(2, 2))

    indx = 0
    if full_id is not None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]

    # =========== Energies & Derivatives ===============
    X, Y = q.get(0, indx), q.get(1, indx)

    H_11 = 0.5 * Kx * (X - X1)**2 + 0.5 * Ky * Y**2
    H_22 = 0.5 * Kx * (X - X2)**2 + 0.5 * Ky * Y**2 + Delta
    expx = math.exp(-alpha * (X - X3)**2)
    expy = math.exp(-beta * Y**2)
    H_12 = gamma * Y * expx * expy

    obj.ham_dia.set(0, 0, H_11 * (1.0 + 0.0j))
    obj.ham_dia.set(1, 1, H_22 * (1.0 + 0.0j))
    obj.ham_dia.set(0, 1, H_12 * (1.0 + 0.0j))
    obj.ham_dia.set(1, 0, H_12 * (1.0 + 0.0j))

    # dH/dX
    dH_11 = Kx * (X - X1)
    dH_22 = Kx * (X - X2)
    dexpx = -2.0 * alpha * (X - X3) * expx
    dH_12 = gamma * Y * dexpx * expy

    obj.d1ham_dia[0].set(0, 0, dH_11 * (1.0 + 0.0j))
    obj.d1ham_dia[0].set(1, 1, dH_22 * (1.0 + 0.0j))
    obj.d1ham_dia[0].set(0, 1, dH_12 * (1.0 + 0.0j))
    obj.d1ham_dia[0].set(1, 0, dH_12 * (1.0 + 0.0j))

    # dH/dY
    dH_11 = Ky * Y
    dH_22 = Ky * Y
    dexpy = -2.0 * beta * Y * expy
    dH_12 = gamma * Y * expx * dexpy + gamma * expx * expy

    obj.d1ham_dia[1].set(0, 0, dH_11 * (1.0 + 0.0j))
    obj.d1ham_dia[1].set(1, 1, dH_22 * (1.0 + 0.0j))
    obj.d1ham_dia[1].set(0, 1, dH_12 * (1.0 + 0.0j))
    obj.d1ham_dia[1].set(1, 0, dH_12 * (1.0 + 0.0j))

    # ========== Second derivatives =================
    # d^2H/dX^2
    d2HdX2_11 = Kx
    d2HdX2_22 = Kx
    # expx = math.exp(-alpha * (X-X3)**2 )
    # dexpx =  -2.0 * alpha * (X-X3) * expx
    d2expx = -2.0 * alpha * (expx + (X - X3) * dexpx)
    d2HdX2_12 = gamma * Y * d2expx * expy  # d^2H11/dX^2

    obj.d2ham_dia[0].set(0, 0, d2HdX2_11 * (1.0 + 0.0j))
    obj.d2ham_dia[0].set(1, 1, d2HdX2_22 * (1.0 + 0.0j))
    obj.d2ham_dia[0].set(0, 1, d2HdX2_12 * (1.0 + 0.0j))
    obj.d2ham_dia[0].set(1, 0, d2HdX2_12 * (1.0 + 0.0j))

    # d^2H/dXdY and d^2H/dYdX
    d2HdXdY_11 = 0.0
    d2HdXdY_22 = 0.0
    # dexpx =  -2.0 * alpha * (X-X3) * expx
    # dH_12 = gamma * Y * dexpx * expy  = dH/dX
    d2HdXdY_12 = gamma * dexpx * (expy + Y * dexpy)
    obj.d2ham_dia[1].set(0, 0, d2HdXdY_11 * (1.0 + 0.0j))
    obj.d2ham_dia[1].set(1, 1, d2HdXdY_22 * (1.0 + 0.0j))
    obj.d2ham_dia[1].set(0, 1, d2HdXdY_12 * (1.0 + 0.0j))
    obj.d2ham_dia[1].set(1, 0, d2HdXdY_12 * (1.0 + 0.0j))

    obj.d2ham_dia[2].set(0, 0, d2HdXdY_11 * (1.0 + 0.0j))
    obj.d2ham_dia[2].set(1, 1, d2HdXdY_22 * (1.0 + 0.0j))
    obj.d2ham_dia[2].set(0, 1, d2HdXdY_12 * (1.0 + 0.0j))
    obj.d2ham_dia[2].set(1, 0, d2HdXdY_12 * (1.0 + 0.0j))

    # d^2H/dY^2
    d2HdY2_11 = Ky
    d2HdY2_22 = Ky
    # dexpy = -2.0 * beta * Y * expy
    d2expy = -2.0 * beta * (Y * dexpy + expy)
    # dH_12 = gamma * Y * expx * dexpy + gamma * expx * expy
    d2HdY2_12 = gamma * expx * (dexpy + Y * d2expy + dexpy)
    obj.d2ham_dia[3].set(0, 0, d2HdY2_11 * (1.0 + 0.0j))
    obj.d2ham_dia[3].set(1, 1, d2HdY2_22 * (1.0 + 0.0j))
    obj.d2ham_dia[3].set(0, 1, d2HdY2_12 * (1.0 + 0.0j))
    obj.d2ham_dia[3].set(1, 0, d2HdY2_12 * (1.0 + 0.0j))

    return obj
