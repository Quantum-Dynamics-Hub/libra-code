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
            * obj.d2ham_dia ( list of 1 CMATRIX(4,4) objects ):
                second derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
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
    d2ham_dia = CMATRIXList();  d2ham_dia.append( CMATRIX(n,n) )
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

    for k in [0]:
        #  d2 Hdia / dR_0 2
        for i in range(n):
            d2ham_dia[k].set(i,i, k_n[i]*(1.0+0.0j) )

    for k in [0]:
        for i in range(n):
            for j in range(n):
                if i!=j:
                    d2ham_dia[k].set(i,j, -2.0*alpha[i][j] * V[i][j] * math.exp(-alpha[i][j] * (x-x_nm[i][j])**2 )*(1.0-2.0*alpha[i][j]*(x-x_nm[i][j])**2)*(1.0+0.0j) )


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.d2ham_dia = d2ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def convert_mtx(X, nx, ny):
    """
    Convert a list of size `nstates x nstates` with objects of size nx x ny
    to a 2D array of the MATRIX(nx, ny) objects
    """
    n_sq = len(X)
    n = int(math.sqrt(n_sq))
    if n*n != n_sq:
        print("All the parameters should be input as lists on nstates x nstates objects. Exiting...\n");
        sys.exit(0)    

    x = []
    for i in range(n):
        xi = []
        for j in range(n):
            # Allocate xi
            xi.append( MATRIX(nx, ny))

            # Populate xi
            for k1 in range(nx):
                for k2 in range(ny):
                    xi[j].set(k1, k2, X[i*n+j][k1*ny+k2])
        x.append(xi)
    return x
    
def polynomial(q, A, B, C, T):
    """
    q = MATRIX(ndof, 1)
    A = MATRIX(ndof, ndof)
    B = MATRIX(ndof, 1)
    C = MATRIX(1,1)
    T = MATRIX(ndof, 1)
    """

    dq = q - T
    F = (0.5 * dq.T() * A * dq + B.T() * dq + C).get(0, 0)
    dF = A * dq + B 

    return F, dF


def Holstein_gen(_q, params, full_id=None):
    """
    This is the generalized Holstein n-level Hamiltonian

    H_nm = F(q; P_nm) * exp(-F(q; p_nm) )\

    P_i = (T_i, A_i, B_i, C_i); i = (n,m)  - pre-exp polynomial parameters
    p_i = (t_i, a_i, b_i, c_i); i = (n,m)  - exp polynomial parameters

    T_i or q_i = MATRIX(ndof, 1)
    A_i or a_i = MATRIX(ndof, ndof)
    B_i or b_i = MATRIX(ndof, 1)
    C_i or c_i = scalar

    F(q; P_i) = 0.5 * (q-T_i)^T * A_i (q-T_i) + B_i^T * (q-T_i) + C_i

    Args:
        _q ( MATRIX(ndof, ntraj) ): coordinates of the particle for all trajectories
        params ( dictionary ): model parameters

        each of the lists below is of size  nstates x nstates, where `nstates` is the number of electronic states

            * **params["A"]** ( list of 2D array of doubles of size ndof x ndof ): parameter for pre-exp polynomial [ units: Ha/Bohr^2 ]
            * **params["B"]** ( list of 1D array of doubles of size ndof ): parameter for pre-exp polynomial [ units: Ha/Bohr]
            * **params["C"]** ( list of doubles ): parameter for pre-exp polynomial [ Ha ]
            * **params["T"]**( list of 1D array of doubles of size ndof): parameter for pre-exp polynomial [ units: Bohr ]
            * **params["a"]** ( list of 2D array of doubles of size ndof x ndof ): parameter for exp polynomial [ units: Ha/Bohr^2 ]
            * **params["b"]** ( list of 1D array of doubles of size ndof ): parameter for exp polynomial [ units: Ha/Bohr]
            * **params["c"]** ( list of doubles ): parameter for exp polynomial [ Ha ]
            * **params["t"]**( list of 1D array of doubles of size ndof): parameter for exp polynomial [ units: Bohr ]

    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(nstates,nstates) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(nstates,nstates) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of ndof CMATRIX(nstates,nstates) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.d2ham_dia ( list of ndof x ndof CMATRIX(nstates,nstates)) objects ):
                second derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of ndof CMATRIX(nstates,nstates) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    critical_params = [ ]

    # Here are parameters for a 2-state systems with 1 nuclear dofs:
    # two parabolas, both depending on 1 nuclear DOF
    # constant coupling
    default_params = {}

    k0, k1 = 0.001, 0.001  # force constants of each diabatic surface
    V = 0.001              # diabatic coupling
    E0, E1 = 0.0, 0.01     # diabatic shifts
    alp0, alp1 = 0.0, 0.0
    x0, x1, x01 = 0.0, 1.0, 0.5
    default_params["A"] = [  [k0],   [0.0], 
                             [0.0],  [k1]
                          ]
    default_params["B"] = [  [0.0],   [0.0],
                             [0.0],   [0.0]
                          ]
    default_params["C"] = [  [E0],   [V],
                             [V],    [E1]
                          ]
    default_params["T"] = [  [x0],   [x01],
                             [x01],   [x1]
                          ]

    default_params["a"] = [  [2.0*alp0],   [0.0],
                             [0.0],  [2.0*alp1]
                          ]
    default_params["b"] = [  [0.0],   [0.0],
                             [0.0],   [0.0]
                          ]
    default_params["c"] = [  [0.0],   [0.0],
                             [0.0],    [0.0]
                          ]
    default_params["t"] = [  [x0],   [x01],
                             [x01],   [x1]
                          ]    

    # Get the parameters
    comn.check_input(params, default_params, critical_params)

    A = params["A"]; 
    B = params["B"]
    C = params["C"]
    T = params["T"]
    a = params["a"]
    b = params["b"]
    c = params["c"]
    t = params["t"]

    n = int(math.sqrt(len(A)))
    ndof = len(B[0])

    #============ Convert to matrices =================
    A = convert_mtx(A, ndof, ndof)  
    B = convert_mtx(B, ndof, 1)
    C = convert_mtx(C, 1, 1)
    T = convert_mtx(T, ndof, 1)

    a = convert_mtx(a, ndof, ndof)
    b = convert_mtx(b, ndof, 1)
    c = convert_mtx(c, 1, 1)
    t = convert_mtx(t, ndof, 1)


    Sdia = CMATRIX(n,n); Sdia.identity()
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(n,n) )
    d2ham_dia = CMATRIXList();  d2ham_dia.append( CMATRIX(n,n) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(n,n) )


    indx = 0
    if full_id!=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]

    q = _q.col(indx)


    obj = tmp()
    #============= Diabatic Hamiltonians ======================
    obj.ham_dia = CMATRIX(n,n)
    obj.d1ham_dia = CMATRIXList()
    for i in range(ndof):
        obj.d1ham_dia.append( CMATRIX(n,n) )

    for i in range(n):
        for j in range(n):
            P1, dP1 = polynomial(q, A[i][j], B[i][j], C[i][j], T[i][j])
            P2, dP2 = polynomial(q, a[i][j], b[i][j], c[i][j], t[i][j])
            h_ij = P1 * math.exp(-P2)
            obj.ham_dia.set(i,j,  h_ij * (1.0+0.0j) )
  
            der = dP1 * math.exp(-P2) - dP2 * h_ij; # MATRIX(ndof, 1)
            for k in range(ndof):
                obj.d1ham_dia[k].set(i, j,  der.get(k, 0) * (1.0+0.0j) )

    #============= Diabatic overlaps ==========================
    obj.ovlp_dia = CMATRIX(n,n); obj.ovlp_dia.identity()

    #============= Diabatic derivative couplings ==============
    obj.dc1_dia = CMATRIXList();
    for i in range(ndof):
        obj.dc1_dia.append( CMATRIX(n,n) )

    #============== Second derivatives of diabatic Ham =========
    obj.d2ham_dia = CMATRIXList();
    for i in range(ndof*ndof):
        obj.d2ham_dia.append( CMATRIX(n,n) )



    return obj
