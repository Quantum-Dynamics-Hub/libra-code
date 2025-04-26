# *********************************************************************************
# * Copyright (C) 2024 Daeho Han and Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# ***********************************************************************************
"""
.. module:: GLVC
   :platform: Unix, Windows
   :synopsis: This module implements the generalized Linear Vibronic Coupling (GLVC) Hamiltonian, which includes
              spin-boson, FMO-type of models, and LVC models
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


def gen_bath_params(_params):
    """
    Generates the parameters (oscillator frequencies and coupling strengths) according to selected
    spectral density

    Args:
        params ( dictionary ): model parameters, can contain:

            * **params["nstates"]** ( int ): the number of electronic states [default: 2 ]
            * **params["num_osc"]** ( int ): the number of oscillators on each of the energy levels [ default: 10 ]

            * **params["spectral_density"]** ( int ) : the type of the bath spectral density to consider:
              - 0: user-specified [default: 0 ]
              - 1: Debye density:  J(w) = (lambda/2)  omega * omega_c/(omega**2 + omega_c**2) [ default: 1]
              - 2: Ohmic density:  J(w) = pi *hbar / 2 * ksi * omega * exp(-omega/omega_c), ksi - Kondo parameter

            * **params["Omega"]** ( double ): the bath characteristic frequency, used with spectral_density == 0 [ units: Ha, default: 0.001 ]
            * **params["lambda"]** ( double ): bath reorganization energy,  used with spectral_density == 0 [ units: Ha, default: 0.001 ]
            * **params["Kondo_ksi"]** (double): Kondo parameter [ default: 1.0 ]

    Notes:
        * the length of the "omega" and "coupl" parameters should be equal to num_osc

    Returns:
        tuple: (list, list)

            * omega ( list of `num_osc` doubles ): bath frequencies [ units: Ha/Bohr^2 ]
            * coupl ( list of `num_osc` doubles ): coupling strengths [ units: Ha/Bohr ]
            * mass ( list of `num_osc` doubles): masses of oscillators [ units: a.u. ]
    """

    params = dict(_params)

    critical_params = []
    default_params = {
        "nstates": 2,
        "num_osc": 2,
        "spectral_density": 1,
        "Omega": 0.001,
        "lambda": 0.001,
        "Kondo_ksi": 1.0}
    comn.check_input(params, default_params, critical_params)

    nstates = params["nstates"]
    num_osc = params["num_osc"]
    spectral_density = params["spectral_density"]
    Omega = params["Omega"]
    Lambda = params["lambda"]
    ksi = params["Kondo_ksi"]

    omega, coupl, mass = [], [], []

    # Debye and Ohmic bath discretizations are according to:
    # Wu, D.; Hu, Z.; Li, J.; Sun, X. Forecasting Nonadiabatic Dynamics Using Hybrid Convolutional
    # Neural Network/Long Short-Term Memory Network. The Journal of Chemical
    # Physics 2021, 155 (22), 224104. https://doi.org/10.1063/5.0073689.

    if spectral_density == 0:  # user-specified frequencies
        pass
    elif spectral_density == 1:  # Debye spectral density
        pref = -1.0 * math.sqrt(2.0 * Lambda / (num_osc + 1.0))
#        omega = [ Omega * math.tan( (k + 0.5)/(2*num_osc) * math.pi ) for k in range(num_osc) ] # frequencies
        omega = [Omega * math.tan(0.5 * math.pi * (1.0 - (k / (num_osc + 1))))
                 for k in range(1, num_osc + 1)]  # frequencies
#        coupl = [ omega[k] * pref for k in range(num_osc) ] # coupling strengths
        coupl = [omega[k] * pref for k in range(num_osc)]  # coupling strengths
        mass = [1.0 for k in range(num_osc)]  # masses
    elif spectral_density == 2:  # Ohmic density
        pref = -1.0 * math.sqrt(Omega * ksi / (num_osc + 1.0))
        omega = [-Omega * math.log(1.0 - (k / (num_osc + 1))) for k in range(1, num_osc + 1)]  # frequencies
        coupl = [omega[k] * pref for k in range(num_osc)]  # coupling strengths
        mass = [1.0 for k in range(num_osc)]  # masses

    else:
        print("Only spectral_density = 0, 1, or 2 are allowed for now. Exiting now...")
        sys.exit(0)

    return omega, coupl, mass


def GLVC(q, _params, full_id=None):
    """
    Generalized Linear Vibronic Coupling Hamiltonian, N-levels, F-dim. (same nuclear DOFs are coupled to all states)

    H = H_s + H_b + H_{sb,1}

    H_s = N x N matrix of numbers - diabatic energies and couplings

    H_b  = sum_{n=0}^{N-1} sum_{f=0}^{F-1} |n> (p_f^2/(2*m_f) + 1/2 * m_f * omega_{n,f}^2 * q_f^2 ) <n|  (the kinetic energy is not included here)
                                                                                 frequencies may in general be different for all states

    H_{sb,1} = sum_{n=0}^{N-1} sum_{f=0}^{F-1} |n> coupl_{n,f} * q_{f}^2 <n|   - diagonal linear coupling terms, couplings may in general
                                                                                 be different for all states

    To be added:

    linear off-diagonal terms

    bilinear diagonal terms

    bilinear off-diagonal terms


    Args:
        q ( MATRIX(ndof, 1) ): coordinates of all the particles, ndof should be taken as N * F, which
        corresponds to having F 1-dimensional oscillators on each of the N levels; organized as in N blocks of F elements

        params ( dictionary ): model parameters, should contain:

            * **params["nstates"]** ( int ): the number of electronic states [default: 2 ]

            * **params["num_osc"]** ( int ): the number of oscillators on each of the energy levels [ units: 2 ]

            * **params["Ham"]**   ( list of lists of double ): H_s values [ default:  [[0.00]*2]*2, units: Ha]

            * **params["coupling_scaling"]** (list of N doubles): linear coupling scaling paramters for each state

            * **params["omega"]** (list of N lists of F doubles): frequencies for all states

            * **params["coupl"]** (list of N lists of F doubles): diagonal linear couplings for all states

            * **params["mass"]** (list of F doubles): masses of nuclear DOFS [ default: [1.0, 1.0], units: a.u.]

    Notes:
        * lenth of the "omega" and "coup" parameters should be equal to num_osc, but the q should have nstates x num_osc rows (dofs)


    Returns:
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(nstates,nstates) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(nstates,nstates) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of ndof CMATRIX(nstats,nstates) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of ndof CMATRIX(nstates,nstates) objects ): derivative coupling in the diabatic basis [ zero ]

    """

    params = dict(_params)

    # Define potential specific constants
    critical_params = ["omega", "coupl"]
    default_params = {"Ham": [[0.00, 0.00], [0.00, 0.00]], "nstates": 2,
                      "num_osc": 2, "coupling_scaling": [1.0, -1.0], "mass": [1.0, 1.0]}
    comn.check_input(params, default_params, critical_params)

    w = params["omega"]
    coupl = params["coupl"]
    nstates = params["nstates"]
    num_osc = params["num_osc"]
    Ham = params["Ham"]
    scl = params["coupling_scaling"]
    mass = params["mass"]

    ndof = q.num_of_rows  # the number of nuclear DOFs

    if (ndof == num_osc):
        pass
    else:
        print(F"The coordinates input should have {ndof} = {num_osc} rows\nExiting now...\n")
        sys.exit(0)

    obj = tmp()
    obj.ham_dia = CMATRIX(nstates, nstates)
    obj.ovlp_dia = CMATRIX(nstates, nstates)
    obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in range(ndof):
        obj.d1ham_dia.append(CMATRIX(nstates, nstates))
        obj.dc1_dia.append(CMATRIX(nstates, nstates))

    indx = 0
    if full_id is not None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]

    # =========== Energies & Derivatives ===============
    for i in range(nstates):
        for j in range(nstates):
            obj.ham_dia.set(i, j, Ham[i][j] * (1.0 + 0.0j))

    for n in range(nstates):
        x = 0.0

        for f in range(num_osc):
            q_f = q.get(f, indx)

            # energy
            w2 = mass[f] * w[n][f]**2

            x = x + 0.5 * w2 * q_f**2 + coupl[n][f] * q_f * scl[n]
            y = w2 * q_f + coupl[n][f] * scl[n]

            # derivative w.r.t. q_nf:
            obj.d1ham_dia[f].add(n, n, y * (1.0 + 0.0j))

        obj.ham_dia.add(n, n, x * (1.0 + 0.0j))

    return obj


def get_GLVC_set1():
    """
    Runeson, J. E.; Manolopoulos, D. E. J. Chem. Phys. 159, 094115, 2023
    FMO, 3-states model

    Args:
        None

    Returns:
        dictionary: params, will contain the parameters:

        * **nstates** (int): the number of electronic states

        * **num_osc** (int): the number of oscillators per electronic state

        * **spectral_density** (int): type of spectral density calculation - Debye

        * **Omega** (double): characteristic frequency of bath [units: Ha]

        * **lambda** (double): bath reorganization energy [units: Ha]

        * **omega** (list of doubles): frequencies of the discretized spectral density [ units: Ha ]

        * **coupl** (list of doubles): electron-nuclear couplings, same for all states, different for each oscillator [ units: Ha ]

        * **mass** (list of doubles): masses of nuclear DOFs [ units: a.u. ]

        * **coupling_scaling** (list of `nstates` doubles): the scaling applied to the overall linear coupling terms for each state [ unitS: N/A ]

        * **beta** (double): inverse temperature, 1/(kB*T) [ units: Ha^-1 ]

        * **Ham** (list of `nstates` lists of `nstates` doubles): diabatic state energies and couplings [ units: Ha ]

    """

    params = {}
    params["nstates"] = 3
    params["num_osc"] = 60
    params["spectral_density"] = 1
    params["Omega"] = 106.14 * units.inv_cm2Ha
    params["lambda"] = 35.0 * units.inv_cm2Ha
    Om, Co, mass = gen_bath_params(params)
    params["omega"] = [list(Om), list(Om), list(Om)]
    params["coupl"] = [list(Co), list(Co), list(Co)]
    params["mass"] = list(mass)
    params["coupling_scaling"] = [1.0, 1.0, 1.0]

    s = units.inv_cm2Ha
    e0 = 12410.0 * s
    e1 = 12530.0 * s
    e2 = 12210.0 * s
    v01 = -87.7 * s
    v02 = 5.5 * s
    v12 = 30.8 * s

    params["Ham"] = [[e0, v01, v02],
                     [v01, e1, v12],
                     [v02, v12, e2]]

    return params


def get_GLVC_set2(indx):
    """
    2-state spin-boson model

    Tempelaar, R.; Reichman, D. R. Generalization of Fewest-Switches Surface Hopping for Coherences.
    The Journal of Chemical Physics 2018, 148 (10), 102309. https://doi.org/10.1063/1.5000843.

    Args:
        indx (int): index of the parameters set:

        - 0 : Figure 2a
        - 1 : Figure 2b
        - 2 : Figure 2c
        - 3 : Figure 2d
        - 4 : Figure 3a
        - 5 : Figure 3b
        - 6 : Figure 3c
        - 7 : Figure 4a
        - 8 : Figure 4b
        - 9 : Figure 5a
        - 10: Figure 5b

    Returns:
        dictionary: params, will contain the parameters:

        * **nstates** (int): the number of electronic states

        * **num_osc** (int): the number of oscillators per electronic state

        * **spectral_density** (int): type of spectral density calculation - Debye

        * **Omega** (double): characteristic frequency of bath [units: Ha]

        * **lambda** (double): bath reorganization energy [units: Ha]

        * **omega** (list of doubles): frequencies of the discretized spectral density [ units: Ha ]

        * **coupl** (list of doubles): electron-nuclear couplings, same for all states, different for each oscillator [ units: Ha ]

        * **mass** (list of doubles): masses of nuclear DOFs [ units: a.u. ]

        * **coupling_scaling** (list of `nstates` doubles): the scaling applied to the overall linear coupling terms for each state [ unitS: N/A ]

        * **beta** (double): inverse temperature, 1/(kB*T) [ units: Ha^-1 ]

        * **Ham** (list of `nstates` lists of `nstates` doubles): diabatic state energies and couplings [ units: Ha ]


    """

    s = 208.5 * units.inv_cm2Ha  # thermal energy at 300 K

    params = {}
    params["nstates"] = 2
    params["num_osc"] = 100
    params["spectral_density"] = 1

    E, V, L, Om, T = 0.0, 0.0, 0.0, 0.0, 0.0

    # Figure 2 - vary lambda
    if indx == 0:
        E, V, L, Om, T = 0.5, 0.5, 0.02, 0.1, 1.0
    elif indx == 1:
        E, V, L, Om, T = 0.5, 0.5, 0.1, 0.1, 1.0
    elif indx == 2:
        E, V, L, Om, T = 0.5, 0.5, 1.0, 0.1, 1.0
    elif indx == 3:
        E, V, L, Om, T = 0.5, 0.5, 5.0, 0.1, 1.0

    # Figure 3 - vary V
    elif indx == 4:
        E, V, L, Om, T = 0.5, 1.0, 1.0, 0.1, 1.0
    elif indx == 5:
        E, V, L, Om, T = 0.5, 0.5, 1.0, 0.1, 1.0
    elif indx == 6:
        E, V, L, Om, T = 0.5, 0.1, 1.0, 0.1, 1.0

    # Figure 4 - vary T, but preserve T * Om
    elif indx == 7:
        E, V, L, Om, T = 0.5, 0.5, 1.0, 0.1, 1.0
    elif indx == 8:
        E, V, L, Om, T = 0.5, 0.5, 1.0, 1.0, 0.1

    # Figure 5 - vary E
    elif indx == 9:
        E, V, L, Om, T = 0.5, 0.5, 1.0, 0.1, 1.0
    elif indx == 10:
        E, V, L, Om, T = 0.0, 0.5, 1.0, 0.1, 1.0

    E = E * s
    V = V * s
    params["Omega"] = Om * s
    params["lambda"] = L * s
    params["beta"] = 1.0 / (T * s)
    OM, CO, mass = gen_bath_params(params)
    params["omega"] = [list(OM), list(OM)]
    params["coupl"] = [list(CO), list(CO)]
    params["mass"] = list(mass)

    params["Ham"] = [[E, V], [V, -E]]
    params["coupling_scaling"] = [1.0, -1.0]

    return params


def get_GLVC_set3():
    """
    5-state, 5 nuclear DOFS chain of coupled oscillators

    General:
    Wang, L.; Prezhdo, O. V. A Simple Solution to the Trivial Crossing Problem in Surface Hopping.
    J. Phys. Chem. Lett. 2014, 5 (4), 713–719. https://doi.org/10.1021/jz500025c.

    with N = 5
    Jain, A.; Alguire, E.; Subotnik, J. E. An Efficient, Augmented Surface Hopping Algorithm That
    Includes Decoherence for Use in Large-Scale Simulations. J. Chem. Theory Comput. 2016, 12 (11),
    5256–5268. https://doi.org/10.1021/acs.jctc.6b00673.

    Args:

    Returns:
        dictionary: params, will contain the parameters:

        * **nstates** (int): the number of electronic states

        * **num_osc** (int): the number of oscillators per electronic state

        * **spectral_density** (int): type of spectral density calculation - Debye

        * **Omega** (double): characteristic frequency of bath [units: Ha]

        * **lambda** (double): bath reorganization energy [units: Ha]

        * **omega** (list of doubles): frequencies of the discretized spectral density [ units: Ha ]

        * **coupl** (list of doubles): electron-nuclear couplings, same for all states, different for each oscillator [ units: Ha ]

        * **mass** (list of doubles): masses of nuclear DOFs [ units: a.u. ]

        * **coupling_scaling** (list of `nstates` doubles): the scaling applied to the overall linear coupling terms for each state [ unitS: N/A ]

        * **beta** (double): inverse temperature, 1/(kB*T) [ units: Ha^-1 ]

        * **Ham** (list of `nstates` lists of `nstates` doubles): diabatic state energies and couplings [ units: Ha ]

    """

    params = {}
    params["nstates"] = 5
    params["num_osc"] = 5
    params["spectral_density"] = 0

    Om = 200.0 * units.inv_cm2Ha
    OM = [Om, Om, Om, Om, Om]
    V = 50.0 * units.inv_cm2Ha
    g = 3091.8 * units.inv_cm2Ha

    # mass should be 1836

    params["Omega"] = Om
    params["lambda"] = 0.0  # doesn't matter
    params["beta"] = 1.0 / (575.5 * units.kB)

    params["omega"] = [list(OM), list(OM), list(OM), list(OM), list(OM)]
    params["coupl"] = [[g, 0.0, 0.0, 0.0, 0.0],
                       [0.0, g, 0.0, 0.0, 0.0],
                       [0.0, 0.0, g, 0.0, 0.0],
                       [0.0, 0.0, 0.0, g, 0.0],
                       [0.0, 0.0, 0.0, 0.0, g]]

    params["Ham"] = [[0.0, V, 0.0, 0.0, 0.0],
                     [V, 0.0, V, 0.0, 0.0],
                     [0.0, V, 0.0, V, 0.0],
                     [0.0, 0.0, V, 0.0, V],
                     [0.0, 0.0, 0.0, V, 0.0]]

    params["mass"] = [1836.0 for i in range(5)]
    params["coupling_scaling"] = [1.0, 1.0, 1.0, 1.0, 1.0]

    return params


def get_GLVC_set4():
    """
    3-state, 3 nuclear DOFS chain of coupled oscillators, top row of Figure 4

    By:
    Wang, L.; Prezhdo, O. V. A Simple Solution to the Trivial Crossing Problem in Surface Hopping.
    J. Phys. Chem. Lett. 2014, 5 (4), 713–719. https://doi.org/10.1021/jz500025c.

    Args:

    Returns:
        dictionary: params, will contain the parameters:

        * **nstates** (int): the number of electronic states

        * **num_osc** (int): the number of oscillators per electronic state

        * **spectral_density** (int): type of spectral density calculation - Debye

        * **Omega** (double): characteristic frequency of bath [units: Ha]

        * **lambda** (double): bath reorganization energy [units: Ha]

        * **omega** (list of doubles): frequencies of the discretized spectral density [ units: Ha ]

        * **coupl** (list of doubles): electron-nuclear couplings, same for all states, different for each oscillator [ units: Ha ]

        * **mass** (list of doubles): masses of nuclear DOFs [ units: a.u. ]

        * **coupling_scaling** (list of `nstates` doubles): the scaling applied to the overall linear coupling terms for each state [ unitS: N/A ]

        * **beta** (double): inverse temperature, 1/(kB*T) [ units: Ha^-1 ]

        * **Ham** (list of `nstates` lists of `nstates` doubles): diabatic state energies and couplings [ units: Ha ]

    """

    params = {}
    params["nstates"] = 3
    params["num_osc"] = 3
    params["spectral_density"] = 0

    Om = 0.01 * units.ev2Ha
    J = 0.005 * units.ev2Ha
    L = 0.1 * units.ev2Ha
    M = 100.0
    K = M * Om**2
    g = math.sqrt(L * K)  # in the paper is called alpha

    params["Omega"] = Om
    OM = [Om, Om, Om]
    params["lambda"] = L
    params["beta"] = 1.0 / (300.0 * units.kB)

    params["omega"] = [list(OM), list(OM), list(OM)]
    params["coupl"] = [[g, 0.0, 0.0],
                       [0.0, g, 0.0],
                       [0.0, 0.0, g]]

    params["Ham"] = [[0.0, J, 0.0],
                     [J, 0.0, J],
                     [0.0, J, 0.0]]

    params["mass"] = [M for i in range(3)]
    params["coupling_scaling"] = [1.0, 1.0, 1.0]

    return params


def get_GLVC_set5(indx):
    """
    3-state spin-boson model

    Bondarenko, A. S.; Tempelaar, R. Overcoming Positivity Violations for Density Matrices in Surface Hopping.
    The Journal of Chemical Physics 2023, 158 (5), 054117. https://doi.org/10.1063/5.0135456.

    Args:
        indx (int): index of the parameters set:


    Returns:
        dictionary: params, will contain the parameters:

        * **nstates** (int): the number of electronic states

        * **num_osc** (int): the number of oscillators per electronic state

        * **spectral_density** (int): type of spectral density calculation - Debye

        * **Omega** (double): characteristic frequency of bath [units: Ha]

        * **lambda** (double): bath reorganization energy [units: Ha]

        * **omega** (list of doubles): frequencies of the discretized spectral density [ units: Ha ]

        * **coupl** (list of doubles): electron-nuclear couplings, same for all states, different for each oscillator [ units: Ha ]

        * **mass** (list of doubles): masses of nuclear DOFs [ units: a.u. ]

        * **coupling_scaling** (list of `nstates` doubles): the scaling applied to the overall linear coupling terms for each state [ unitS: N/A ]

        * **beta** (double): inverse temperature, 1/(kB*T) [ units: Ha^-1 ]

        * **Ham** (list of `nstates` lists of `nstates` doubles): diabatic state energies and couplings [ units: Ha ]


    """

    s = 208.5 * units.inv_cm2Ha  # thermal energy at 300 K

    params = {}
    params["nstates"] = 3
    params["num_osc"] = 100  # not described in the paper, but let's assume it is the same as in their 2-state model
    params["spectral_density"] = 1

    E1, E2, E3, V, L, Om, T = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if indx == 0:  # Figure 1 - homogeneous trimer, adiabatic regime at low T
        E1, E2, E3, V, L, Om, T = 0.0, 0.0, 0.0, 1.0, 0.005, 1.0, 0.1
    elif indx == 1:  # Figre 2 - biased trimer
        E1, E2, E3, V, L, Om, T = 0.0, -0.5, -1.0, 1.0, 0.005, 0.1, 1.0
    elif indx == 2:  # Figure 6 - homogeneous trimer, adiabatic regime at low T, larger reorg. energy
        E1, E2, E3, V, L, Om, T = 0.0, 0.0, 0.0, 1.0, 0.025, 1.0, 0.1
    elif indx == 3:  # top figure first column on page S4 of SI
        E1, E2, E3, V, L, Om, T = 0.0, 0.0, 0.0, 1.0, 0.005, 0.1, 1.0
    elif indx == 4:  # top figure third column on page S4 of SI
        E1, E2, E3, V, L, Om, T = 0.0, 0.0, 0.0, 1.0, 0.025, 0.1, 1.0  # temporary change
    elif indx == 5:  # bottom figure first column on page S4 of SI
        E1, E2, E3, V, L, Om, T = 0.0, 0.0, 0.0, 0.1, 0.005, 0.1, 1.0
    elif indx == 6:  # bottom figure third column on page S4 of SI
        E1, E2, E3, V, L, Om, T = 0.0, 0.0, 0.0, 0.1, 0.25, 0.1, 1.0

    E1 = E1 * s
    E2 = E2 * s
    E3 = E3 * s
    V = V * s
    params["Omega"] = Om * s
    params["lambda"] = L * s
    params["beta"] = 1.0 / (T * s)
    OM, CO, mass = gen_bath_params(params)
    params["omega"] = [list(OM), list(OM), list(OM)]
    params["coupl"] = [list(CO), list(CO), list(CO)]
    params["mass"] = list(mass)

    params["Ham"] = [[E1, V, 0.0],
                     [V, E2, V],
                     [0.0, V, E3]]
    params["coupling_scaling"] = [1.0, 1.0, 1.0]

    return params


def get_GLVC_set6(indx, scl=1.0):
    """
    2-state spin-boson model

    Mannouch, J. R.; Richardson, J. O. A Mapping Approach to Surface Hopping. The Journal of Chemical Physics 2023,
    158 (10), 104111. https://doi.org/10.1063/5.0139734.

    Args:
        indx (int): index of the parameters set:

        - 0 : Figure 10a
        - 1 : Figure 10b
        - 2 : Figure 10c

        scl (double): model scaling, value of 1 is chosen so that
           t * \\Delta = 20 corresponds to t = 2000 a.u. that is \\Delta = 0.01 a.u.


    Returns:
        dictionary: params, will contain the parameters:

        * **nstates** (int): the number of electronic states

        * **num_osc** (int): the number of oscillators per electronic state

        * **spectral_density** (int): type of spectral density calculation - Debye

        * **Omega** (double): characteristic frequency of bath [units: Ha]

        * **lambda** (double): bath reorganization energy [units: Ha]

        * **omega** (list of doubles): frequencies of the discretized spectral density [ units: Ha ]

        * **coupl** (list of doubles): electron-nuclear couplings, same for all states, different for each oscillator [ units: Ha ]

        * **mass** (list of doubles): masses of nuclear DOFs [ units: a.u. ]

        * **coupling_scaling** (list of `nstates` doubles): the scaling applied to the overall linear coupling terms for each state [ unitS: N/A ]

        * **beta** (double): inverse temperature, 1/(kB*T) [ units: Ha^-1 ]

        * **Ham** (list of `nstates` lists of `nstates` doubles): diabatic state energies and couplings [ units: Ha ]


    """

    s = scl  # scaling of the energy units
    delta = 0.0005 * s

    params = {}
    params["nstates"] = 2
    params["num_osc"] = 100
    params["spectral_density"] = 1

    E = delta
    L = 0.5 * delta
    V = delta
    beta = 1.0
    Om = 0.0
    # Figure 10 - vary temperature and Omega
    if indx == 0:
        Om, beta = 0.2 * delta, 0.25 / delta
    elif indx == 1:
        Om, beta = delta, 0.25 / delta
    elif indx == 2:
        Om, beta = delta, 5.0 / delta

    params["Omega"] = Om
    params["lambda"] = L
    params["beta"] = beta
    OM, CO, mass = gen_bath_params(params)
    params["omega"] = [list(OM), list(OM)]
    params["coupl"] = [list(CO), list(CO)]
    params["mass"] = list(mass)

    params["Ham"] = [[E, V], [V, -E]]
    params["coupling_scaling"] = [1.0, -1.0]

    return params
