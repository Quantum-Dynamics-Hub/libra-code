#*********************************************************************************
#* Copyright (C) 2018-2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: fgr_py
   :platform: Unix, Windows
   :synopsis: 
       This module implements the auxiliary Python functions for dealing with 
       various kinds of FGR calculations (as implemented in C++)

.. moduleauthor:: Alexey V. Akimov

"""


import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import common_utils as comn
import units

def run_NEFGRL_populations(omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, params):
    """

    Noneq FGR in non-Condon case (linear coupling) using normal modes
     k(t') = 2/(hbar^2) * Re { int_0^t' { dtau * C(tau) }  }
    
    Compute the non-equilibrium FGR populations as a function of time

    Args:
        omega_DA ( double ): energy gap between donor and acceptor states E(donor) - E(acceptor) [units: a.u.]
        V ( double ): electronic coupling between the two states [units: a.u.]
        omega_nm ( list of doubles ): frequencies of the bath normal modes [units: a.u.]
        gamma_nm ( list of doubles ): couplings of the bath modes to the quantum system/primary mode [units: a.u.]
        req_nm ( list of doubles ): displacements of the normal modes from their equilibrium positions upon the
            charge transfer 
        shift_NE ( list of doubles ): non-equilibrium displacements of the normal modes

        params ( dictionary ): parameters controlling the execution of the calculations

            * **params["tmax"]** ( double ): time since the initial photoexcitation [units: a.u. of time]

            * **params["dt"]** ( double ): integration timestep for the forward propagation [units: a.u. of time]

            * **params["dtau"]** ( double ): integration timestep for backward propagation, 
                used to take the integral above (tau = dtau*n < t, where n is integer)

            * **params["method"]** ( int ): flag that specifies which method to use:

                - 0: Exact quantum result
                - 1: LSC
                - 2: CAV
                - 3: CD
                - 4: W0
                - 5: Marcus

            * **params["dyn_type"]** ( int ): flag that selects Condon (0) vs. non-Condon(1) approximation

            * **params["Temperature"]** ( double ): temperature of the bath [ units: K ]


  Returns: the matrix with: current time, instantaneous rate, population on the donor state

        *     
    """


    critical_params = [  ] 
    default_params = {  "method":0, "dyn_type":0, 
                        "Temperature":300.0, 
                        "dtau":1.0, "tmax":10.0, "dt":1.0, 
                        "do_output":False, "filename":"FGR.txt"
                     }
    comn.check_input(params, default_params, critical_params)


    T = params["Temperature"]

    
    f = open(filename, "w")
    f.close()

    nsteps = int(tmax/dt)+1
    summ, P, k = 0.0, 1.0, 0.0  # probability of donor state


    for step in xrange(nsteps):

        t = step*dt

        f = open(filename, "a")
        f.write("%8.5f  %8.5f  %8.5f \n" % (t, k, P))
        f.close()
 
        # k = k(t')
        k = NEFGRL_rate(t, omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, method, beta, dyn_type, dtau)
        
        summ += k * dt; 
        P = math.exp(-summ)  # exp(- int dt' k(t'))



def run_Test1(omega_DA, V, omega_nm, gamma_nm, req_nm, shift_NE, method, beta, dyn_type, t, dtau, filename):


    _tau = []
    _argg_re, _argg_im = [], []
    _lin_re, _lin_im = [], []
    _C_re, _C_im = [], []
    _int_re, _int_im = [], []

    f = open(filename, "w")
    f.close()

    nsteps = int(t/dtau)
    nomega = len(omega_nm)

    """
     int_0_t { C(t,tau) dtau }  for a fixed t
    """

    integ = 0.0+0.0j

    for step in xrange(nsteps):
        tau = step*dtau
        argg, lin = 0.0+0.0j, 0.0+0.0j 

        for w in range(nomega):
            if method==0:
                argg = argg + Integrand_NE_exact(t, tau, omega_DA, omega_nm[w], shift_NE[w], req_nm[w], beta)
                lin = lin + Linear_NE_exact(t, tau, gamma_nm[w], omega_nm[w], shift_NE[w], req_nm[w], beta)
            elif method==1:
                argg = argg + Integrand_NE_LSC(t, tau, omega_DA, omega_nm[w], shift_NE[w], req_nm[w], beta)
                lin = lin + Linear_NE_LSC(t, tau, gamma_nm[w], omega_nm[w], shift_NE[w], req_nm[w], beta)
            elif method==2:
                argg = argg + Integrand_NE_CAV(t, tau, omega_DA, omega_nm[w], shift_NE[w], req_nm[w], beta)
                lin = lin + Linear_NE_CAV(t, tau, gamma_nm[w], omega_nm[w], shift_NE[w], req_nm[w], beta)
            elif method==3:
                argg = argg + Integrand_NE_CD(t, tau, omega_DA, omega_nm[w], shift_NE[w], req_nm[w], beta)
                lin = lin + Linear_NE_CD(t, tau, gamma_nm[w], omega_nm[w], shift_NE[w], req_nm[w], beta)
            elif method==4:
                argg = argg + Integrand_NE_W0(t, tau, omega_DA, omega_nm[w], shift_NE[w], req_nm[w], beta)
                lin = lin + Linear_NE_W0(t, tau, gamma_nm[w], omega_nm[w], shift_NE[w], req_nm[w], beta)
            elif method==5:
                argg = argg + Integrand_NE_Marcus(t, tau, omega_DA, omega_nm[w], shift_NE[w], req_nm[w], beta)
                lin = lin + Linear_NE_Marcus(t, tau, gamma_nm[w], omega_nm[w], shift_NE[w], req_nm[w], beta)

        C = 0.0+0.0j
        if dyn_type==0:
            C = cmath.exp(argg) * V * V
        elif dyn_type==1:
            C = cmath.exp(argg) * lin

        integ = integ + C*dtau

        f = open(filename, "a")
        f.write("%8.5f  %8.5f  %8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f  %8.5f\n" 
                % (tau, argg.real, argg.imag, lin.real, lin.imag, C.real, C.imag, integ.real, integ.imag))
        f.close()


        _tau.append(tau)
        _argg_re.append(argg.real)
        _argg_im.append(argg.imag)
        
        _lin_re.append(lin.real)
        _lin_im.append(lin.imag)
        
        _C_re.append(C.real)
        _C_im.append(C.imag)
        
        _int_re.append(integ.real)
        _int_im.append(integ.imag)


        return _tau, _argg_re, _argg_im, _lin_re, _lin_im, _C_re, _C_im, _int_re, _int_im



