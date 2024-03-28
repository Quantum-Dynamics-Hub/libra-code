#*********************************************************************************                     
#* Copyright (C) 2018-2024 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
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
    Generalized Linear Vibronic Coupling Hamiltonian, N-level, F-dim. for each of N levels problem:

    H = H_s + H_b + H_sb

    H_s = N x N matrix of numbers - diabatic energies and couplings

    H_b  = sum_{n=0}^{N-1} sum_{f=0}^{F-1} |n> (1/2 * p_{n,f}^2 + 1/2 * omega_f^2 * q_{n,f}^2 ) <n|  (the kinetic energy is not included here)

    H_sb = sum_{n=0}^{N-1} sum_{f=0}^{F-1} |n> coupl_{f} * q_{n,f}^2 <n| 
    

    Args:
        q ( MATRIX(ndof, 1) ): coordinates of all the particles, ndof should be taken as N * F, which
        corresponds to having F 1-dimensional oscillators on each of the N levels; organized as in N blocks of F elements

        params ( dictionary ): model parameters, should contain:

            * **params["nstates"]** ( int ): the number of electronic states [default: 2 ]
            * **params["num_osc"]** ( int ): the number of oscillators on each of the energy levels [ units: 10 ]

            * **params["Ham"]**   ( list of lists of double ): H_s values [ default:  [[0.00]*2]*2, units: Ha]

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
    default_params = { "Ham": [ [0.00, 0.00], [0.00, 0.00] ], "nstates":2, "num_osc":10 }
    comn.check_input(params, default_params, critical_params)

    w = params["omega"]
    c = params["coupl"]
    nstates = params["nstates"]
    num_osc = params["num_osc"]
    Ham = params["Ham"]

    ndof = q.num_of_rows  # the number of nuclear DOFs  

    if(ndof == nstates * num_osc):
        pass
    else:
        print(F"The coordinates input should have {nstates} x {num_osc} = {nstates * num_osc} rows\nExiting now...\n")
        sys.exit(0)


    obj = tmp()
    obj.ham_dia = CMATRIX(nstates, nstates)
    obj.ovlp_dia = CMATRIX(nstates,nstates);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in range(0,ndof):
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
            q_nf = q.get(n*num_osc+f, indx)

            # energy
            x = x + 0.5 * w[f] * w[f]* q_nf * q_nf + coupl[f] * q_nf

            y = w[f] * q_nf + coupl[f]

            # derivative w.r.t. q_nf:
            obj.d1ham_dia[n*num_osc + f].add(n,n, y*(1.0+0.0j))

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

    """

    params = {}
    params["nstates"] = 3 
    params["num_osc"] = 60 
    params["spectral_density"] = 1
    params["Omega"] = 106.14 * units.inv_cm2Ha
    params["lambda"] = 35.0 * units.inv_cm2Ha
    params["omega"], params["coupl"] = gen_bath_params(params)

    s = units.inv_cm2Ha
    e0 = 12410.0*s
    e1 = 12530.0*s
    e2 = 12210.0*s
    v01 = -87.7*s
    v02 = 5.5*s
    v12 = 30.8

    params["Ham"] =[ [e0, v01, v02], 
                     [v01, e2, v12],
                     [v02, v12, e3] ]
    
    return params


