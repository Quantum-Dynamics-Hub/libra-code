#*********************************************************************************                     
#* Copyright (C) 2024 Daeho Han and Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#***********************************************************************************
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

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
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
              - 2: Ohmic density (not implemented yet)

            * **params["Omega"]** ( double ): the bath characteristic frequency, used with spectral_density == 0 [ units: Ha, default: 0.001 ]
            * **params["lambda"]** ( double ): bath reorganization energy,  used with spectral_density == 0 [ units: Ha, default: 0.001 ]

    Notes: 
        * the length of the "omega" and "coupl" parameters should be equal to num_osc 

    Returns:
        tuple: (list, list)

            * omega ( list of `num_osc` doubles ): bath frequencies [ units: Ha/Bohr^2 ]
            * coupl ( list of `num_osc` doubles ): coupling strengths [ units: Ha/Bohr ]
    """

    params = dict(_params)

    critical_params = [ ]
    default_params = {"nstates":2, "num_osc":10, "spectral_density":0, "Omega":0.001, "lambda":0.001 }
    comn.check_input(params, default_params, critical_params)

    nstates = params["nstates"]
    num_osc = params["num_osc"]
    spectral_density = params["spectral_density"]
    Omega = params["Omega"]
    Lambda = params["lambda"]

    omega, coupl = [], []

    if spectral_density == 0: # user-specified frequencies
        pass
    elif spectral_density == 1: # Debye spectral density
        pref = math.sqrt(2.0*Lambda/num_osc)
        omega = [ Omega * math.tan( (k + 0.5)/(2*num_osc) * math.pi ) for k in range(num_osc) ] # frequencies
        coupl = [ omega[k] * pref for k in range(num_osc) ] # coupling strengths

    else:
        print("Only spectral_density = 0 or 1 are allowed. Exiting now...")
        sys.exit(0)

    return omega, coupl



def GLVC(q, _params, full_id=None):
    """
    Generalized Linear Vibronic Coupling Hamiltonian, N-levels, F-dim. (same nuclear DOFs are coupled to all states)

    H = H_s + H_b + H_{sb,1} 

    H_s = N x N matrix of numbers - diabatic energies and couplings

    H_b  = sum_{n=0}^{N-1} sum_{f=0}^{F-1} |n> (1/2 * p_f^2 + 1/2 * omega_{n,f}^2 * q_f^2 ) <n|  (the kinetic energy is not included here)
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

            * **params["num_osc"]** ( int ): the number of oscillators on each of the energy levels [ units: 10 ]

            * **params["Ham"]**   ( list of lists of double ): H_s values [ default:  [[0.00]*2]*2, units: Ha]

            * **params["coupling_scaling"]** (list of N doubles): linear coupling scaling paramters for each state

            * **params["omega"]** (list of N lists of F doubles): frequencies for all states

            * **params["coupl"]** (list of N lists of F doubles): diagonal linear couplings for all states

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
    critical_params = [ "omega", "coupl" ] 
    default_params = { "Ham": [ [0.00, 0.00], [0.00, 0.00] ], "nstates":2, "num_osc":10, "coupling_scaling":[1.0, -1.0] }
    comn.check_input(params, default_params, critical_params)

    w = params["omega"]
    coupl = params["coupl"]
    nstates = params["nstates"]
    num_osc = params["num_osc"]
    Ham = params["Ham"]
    scl = params["coupling_scaling"]
  
    ndof = q.num_of_rows  # the number of nuclear DOFs  

    if(ndof == num_osc):
        pass
    else:
        print(F"The coordinates input should have {ndof} = {num_osc} rows\nExiting now...\n")
        sys.exit(0)


    obj = tmp()
    obj.ham_dia = CMATRIX(nstates, nstates)
    obj.ovlp_dia = CMATRIX(nstates,nstates);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()


    for i in range(ndof):
        obj.d1ham_dia.append( CMATRIX(nstates,nstates) )
        obj.dc1_dia.append( CMATRIX(nstates,nstates) )

    indx = 0
    if full_id !=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]


    #=========== Energies & Derivatives ===============
    for i in range(nstates):
        for j in range(nstates):
            obj.ham_dia.set(i,j, Ham[i][j]*(1.0+0.0j))

     
    for n in range(nstates):
        x = 0.0
        
        for f in range(num_osc):
            q_f = q.get(f, indx)
        
            # energy
            w2 = w[n][f]**2

            x = x + 0.5 * w2 * q_f**2 + coupl[n][f] * q_f * scl[n]
            y = w2 * q_f + coupl[n][f] * scl[n]

            # derivative w.r.t. q_nf:
            obj.d1ham_dia[f].add(n,n, y*(1.0+0.0j))
        
        obj.ham_dia.add(n,n, x * (1.0+0.0j))
         

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

    """

    params = {}
    params["nstates"] = 3 
    params["num_osc"] = 60 
    params["spectral_density"] = 1
    params["Omega"] = 106.14 * units.inv_cm2Ha
    params["lambda"] = 35.0 * units.inv_cm2Ha
    Om, Co = gen_bath_params(params)
    params["omega"] = [ list(Om), list(Om), list(Om) ]
    params["coupl"] = [ list(Co), list(Co), list(Co) ]
    params["coupling_scaling"] = [1.0, 1.0, 1.0]

    s = units.inv_cm2Ha
    e0 = 12410.0*s
    e1 = 12530.0*s
    e2 = 12210.0*s
    v01 = -87.7*s
    v02 = 5.5*s
    v12 = 30.8*s

    params["Ham"] =[ [e0, v01, v02], 
                     [v01, e1, v12],
                     [v02, v12, e2] ]
    
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

        * **beta** (double): inverse of thermal energy 1/kT [units: Ha^-1]

    """

    s = 208.5 * units.inv_cm2Ha # thermal energy at 300 K

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
    params["beta"] = 1.0/(T*s)
    OM, CO = gen_bath_params(params)
    params["omega"] = [ list(OM), list(OM) ]
    params["coupl"] = [ list(CO), list(CO) ]

    params["Ham"] = [ [E, V], [V, -E] ]
    params["coupling_scaling"] = [1.0, -1.0]

    return params