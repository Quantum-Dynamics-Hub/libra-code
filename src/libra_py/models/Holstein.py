#*********************************************************************************
#* Copyright (C) 2020 Story Temen, Alexey V. Akimov
#* Copyright (C) 2018-2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#***********************************************************************************
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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units


class tmp:
    pass


def Holstein_uncoupled(q, params):
    """
    Implementation of a generic Holstein Hamiltonian.

    Args:
        q ( MATRIX(ndof, 1) ): coordinates of the classical particles, ndof is an
            arbitrary number of degrees of freedom (e.g. 3N, where N is the number of particles)
        params ( dictionary ): the parameters of the Hamiltonian, should contain:

            * **params["k_harmonic"]** ( double ) [ units: Ha/Bohr^2 ]
            * **params["el-phon_coupling"]** ( double ) [ units: Ha/Bohr ]
            * **params["site_coupling"]** ( double ): electronic coupling between nearby sites [ units: Ha ]
            * **params["is_periodic"]** ( Boolean ): whether the first and last sites are connected

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(ndof,ndof) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(ndof,ndof) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of ndof CMATRIX(ndof,ndof) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 2 CMATRIX(ndof,ndof) objects ): derivative coupling in the diabatic basis [ zero ]



    """

    critical_params = [ "k_harmonic", "el-phon_coupling", "site_coupling", "is_periodic" ]
    default_params = { }
    comn.check_input(params, default_params, critical_params)

    k = params["k_harmonic"]  # force constant
    alpha = params["el-phon_coupling"] # local electron-phonon coupling
    V = params["site_coupling"]  # diabatic coupling between the adjacent sites
    is_periodic = params["is_periodic"]



    ndof = q.num_of_rows  # the number of nuclear DOFs
    N = ndof              # in this case, each site has one nuclear DOF

    obj = tmp()
    obj.ham_dia = CMATRIX(N,N)
    obj.ovlp_dia = CMATRIX(N,N);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in range(0,N):
        obj.d1ham_dia.append( CMATRIX(N,N) )
        obj.dc1_dia.append( CMATRIX(N,N) )


    #=========== Energies & Derivatives ===============
    for i in range(0,N-1):
        obj.ham_dia.set(i,i+1, -V)
        obj.ham_dia.set(i+1,i, -V)

    if is_periodic==True:
        obj.ham_dia.set(0,N-1, -V)
        obj.ham_dia.set(N-1,0, -V)


    Epot = 0.0
    for i in range(0,N):
        x = q.get(i);  Epot += x*x
    Epot = 0.5*k*Epot

    for i in range(0,N):
        x = q.get(i);
        obj.ham_dia.set(i,i, Epot + alpha*x*(1.0+0.0j))


    for i in range(0,N):
        for n in range(0,N):
            x = q.get(n);

            if(n==0):
                obj.d1ham_dia[n].set(i,i, (k*x + alpha)*(1.0+0.0j))  #  dH(i,i)/dx_n
            else:
                obj.d1ham_dia[n].set(i,i, k*x*(1.0+0.0j))            #  dH(i,i)/dx_n

    return obj



def get_Holstein_set1():
    """

    Parameters from:
    Qiu, J.; Bai, X.; Wang, L. Crossing Classified and Corrected Fewest Switches Surface Hopping.
    J. Phys. Chem. Lett. 2018, 9, 4319-4325.

    Args:
        None

    Returns:
        dictionary: params, will contain the parameters:

            * **params["k_harmonic"]** ( double ) [ units: Ha/Bohr^2 ]
            * **params["el-phon_coupling"]** ( double ) [ units: Ha/Bohr ]
            * **params["mass"]** ( double ): mass of the particles [ units: a.u. of mass ]
            * **params["site_coupling"]** ( double ): electronic coupling between nearby sites [ units: Ha ]
            * **params["is_periodic"]** ( Boolean ): whether the first and last sites are connected

    """

    params = {}
    params["k_harmonic"]  = 14500.0 * units.amu/(units.ps2au * units.ps2au)
    params["el-phon_coupling"] = 3500.0 * (units.inv_cm2Ha/units.Angst)
    params["mass"] = 250.0 * units.amu
    params["site_coupling"] = 10.0 * units.inv_cm2Ha
    params["is_periodic"] = False

    return params






def Holstein2(q, params, full_id):
    """
    n-state model

    H_nn = E_n + 0.5*k*(x-x_n)^2

    H_n,n+1 = H_n+1,n = V,
    H_n,m = 0, otherwise

    Args:
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
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

    critical_params = ["E_n", "x_n", "k_n" ]
    default_params = { "V":0.001 }
    comn.check_input(params, default_params, critical_params)

    E_n = params["E_n"]
    x_n = params["x_n"]
    k_n = params["k_n"]
    V = params["V"]

    n = len(E_n)

    Hdia = CMATRIX(n,n)
    Sdia = CMATRIX(n,n)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(n,n) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(n,n) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    x = q.col(indx).get(0)

    Sdia.identity()

    for i in range(n):
        Hdia.set(i,i,  (E_n[i] + 0.5*k_n[i]*(x - x_n[i])**2) * (1.0+0.0j) )

    for i in range(n-1):
        Hdia.set(i,i+1,  V * (1.0+0.0j) )
        Hdia.set(i+1,i,  V * (1.0+0.0j) )

    for k in [0]:
        #  d Hdia / dR_0
        for i in range(n):
            d1ham_dia[k].set(i,i, (k_n[i] * (x - x_n[i]))*(1.0+0.0j) )


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def Holstein3(q, params, full_id):
    """
    n-state model

    H_nn = E_n + 0.5*k*(x-x_n)^2

    H_n,n+m = H_n+m,n = V_m


    Args:
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["E_n"]** ( list of doubles ):  [ default: [0.0, 0.001, 0.001, 0.001], units: Ha]
            * **params["x_n"]** ( list of doubles ):  [ default: [0.0, 1.0, 1.0, 1.0], units: Bohr]
            * **params["k_n"]** ( list of doubles ):  [ default: [0.001, 0.001, 0.001, 0.001], units: Ha/Bohr^2]
            * **params["V_n"]** ( list of doubles ):  The coupling between state i and j, where i and j are n+1 states away from each other. ie H_0,1 = V_0; H_0,2 = V_1
            [ default: [0.001, 0.0001, 0] units: Ha]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(4,4) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(4,4) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(4,4) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(4,4) objects ): derivative coupling in the diabatic basis [ zero ]

    Example:

        "H_nn":[0.0, 0.002, 0.002, 0.002]
        "V_n":[0.001, 0.0001, 0]
        Ham:    0.0     0.001   0.0001   0
                        0.002   0.001    0.0001
                                0.002    0.001
                                         0.002
    """

    critical_params = ["E_n", "x_n", "k_n", "V_n" ]
    default_params = {"E_n":[0.0, 0.001, 0.001, 0.001], "x_n":[0.0, 1.0, 1.0, 1.0], "k_n":[0.001, 0.001, 0.001, 0.001], "V_n":[0.001, 0, 0] }
    comn.check_input(params, default_params, critical_params)

    E_n = params["E_n"]
    x_n = params["x_n"]
    k_n = params["k_n"]
    V_n = params["V_n"]

    n = len(E_n)

    Hdia = CMATRIX(n,n)
    Sdia = CMATRIX(n,n)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(n,n) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(n,n) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    x = q.col(indx).get(0)

    Sdia.identity()

    for i in range(n):
        Hdia.set(i,i,  (E_n[i] + 0.5*k_n[i]*(x - x_n[i])**2) * (1.0+0.0j) )

        for k in range(n):	# k = 'distance' between states
            if k != i:
                Hdia.set(i,k, V_n[abs(i-k)-1] * (1.0+0.0j) )

    for k in [0]:
        #  d Hdia / dR_0
        for i in range(n):
            d1ham_dia[k].set(i,i, (k_n[i] * (x - x_n[i]))*(1.0+0.0j) )

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def Holstein4(q, params, full_id):
    """
    n-state model

    H_nn = E_n + 0.5*k*(x-x_n)^2
    H_n,m = V_n,m, 

    Args:
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["E_n"]** ( list of doubles ):  [ default: [0.0, 0.001, 0.001, 0.001], units: Ha]
            * **params["x_n"]** ( list of doubles ):  [ default: [0.0, 1.0, 1.0, 1.0], units: Bohr]
            * **params["k_n"]** ( list of doubles ):  [ default: [0.001, 0.001, 0.001, 0.001], units: Ha/Bohr^2]
            * **params["V"]**   ( list of lists of double ):  [ default:  [[0.001]*4]*4, units: Ha]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(4,4) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(4,4) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(4,4) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(4,4) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = ["E_n", "x_n", "k_n" ]
    default_params = { "V": [ [0.001, 0.001, 0.001, 0.001], [0.001, 0.001, 0.001, 0.001],
                              [0.001, 0.001, 0.001, 0.001], [0.001, 0.001, 0.001, 0.001] ]
                     }
    comn.check_input(params, default_params, critical_params)

    E_n = params["E_n"]
    x_n = params["x_n"]
    k_n = params["k_n"]
    V = params["V"]

    n = len(E_n)

    Hdia = CMATRIX(n,n)
    Sdia = CMATRIX(n,n)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(n,n) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(n,n) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    x = q.col(indx).get(0)

    Sdia.identity()

    for i in range(n):
        Hdia.set(i,i,  (E_n[i] + 0.5*k_n[i]*(x - x_n[i])**2) * (1.0+0.0j) )

    for i in range(n):
        for j in range(n):
            if i!=j:
                Hdia.set(i,j,  V[i][j] * (1.0+0.0j) )

    for k in [0]:
        #  d Hdia / dR_0
        for i in range(n):
            d1ham_dia[k].set(i,i, (k_n[i] * (x - x_n[i]))*(1.0+0.0j) )


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj




def Holstein5(q, params, full_id):
    """
    n-state model

    H_nn = E_n + 0.5*k*(x-x_n)^2
    H_n,m = V_n,m * exp(- alp_n,m * (x-x_nm)^2 )

    Args:
        q ( MATRIX(1,1) ): coordinates of the particle, ndof = 1
        params ( dictionary ): model parameters

            * **params["E_n"]** ( list of doubles ):  [ default: [0.0, 0.001, 0.001, 0.001], units: Ha]
            * **params["x_n"]** ( list of doubles ):  [ default: [0.0, 1.0, 1.0, 1.0], units: Bohr]
            * **params["k_n"]** ( list of doubles ):  [ default: [0.001, 0.001, 0.001, 0.001], units: Ha/Bohr^2]
            * **params["V"]**   ( list of lists of double ):  [ default:  [[0.001]*4]*4, units: Ha]
            * **params["alpha"]** ( list of lists of double ):  [ default:  [[0.0]*4]*4, units: 1/Bohr^2]
            * **params["x_nm"]**  ( list of lists of double ):  [ default:  [[0.0]*4]*4, units: Bohr]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(4,4) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(4,4) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of 1 CMATRIX(4,4) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of 1 CMATRIX(4,4) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = ["E_n", "x_n", "k_n" ]
    default_params = { "V": [ [0.001, 0.001, 0.001, 0.001], [0.001, 0.001, 0.001, 0.001],
                              [0.001, 0.001, 0.001, 0.001], [0.001, 0.001, 0.001, 0.001] ],
                       "alpha": [ [0.00, 0.00, 0.00, 0.00], [0.00, 0.00, 0.00, 0.00],
                                  [0.00, 0.00, 0.00, 0.00], [0.00, 0.00, 0.00, 0.00] ],
                       "x_nm": [ [0.00, 0.00, 0.00, 0.00], [0.00, 0.00, 0.00, 0.00],
                                 [0.00, 0.00, 0.00, 0.00], [0.00, 0.00, 0.00, 0.00] ],
                     }
    comn.check_input(params, default_params, critical_params)

    E_n = params["E_n"]
    x_n = params["x_n"]
    k_n = params["k_n"]
    V = params["V"]
    alpha = params["alpha"]
    x_nm = params["x_nm"]

    n = len(E_n)

    Hdia = CMATRIX(n,n)
    Sdia = CMATRIX(n,n)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(n,n) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(n,n) )

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    x = q.col(indx).get(0)

    Sdia.identity()

    for i in range(n):
        Hdia.set(i,i,  (E_n[i] + 0.5*k_n[i]*(x - x_n[i])**2) * (1.0+0.0j) )

    for i in range(n):
        for j in range(n):
            if i!=j:
                Hdia.set(i,j,  V[i][j] * math.exp(-alpha[i][j] * (x-x_nm[i][j])**2 ) * (1.0+0.0j) )

    for k in [0]:
        #  d Hdia / dR_0
        for i in range(n):
            d1ham_dia[k].set(i,i, (k_n[i] * (x - x_n[i]))*(1.0+0.0j) )

    for k in [0]:
        for i in range(n):
            for j in range(n):
                if i!=j:
                    d1ham_dia[k].set(i,j,  -2.0*alpha[i][j] * (x-x_nm[i][j]) * V[i][j] * math.exp(-alpha[i][j] * (x-x_nm[i][j])**2 ) * (1.0+0.0j) )



    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
