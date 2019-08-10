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
.. module:: models_Holstein
   :platform: Unix, Windows
   :synopsis: This module implements the Henon-Heiles Hamiltonians
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
#import common_utils as comn
import util.libutil as comn
import units


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

    for i in xrange(N):
        obj.d1ham_dia.append( CMATRIX(N,N) )
        obj.dc1_dia.append( CMATRIX(N,N) )


    #=========== Energies & Derivatives ===============
    for i in xrange(N-1):
        obj.ham_dia.set(i,i+1, -V)
        obj.ham_dia.set(i+1,i, -V)

    if is_periodic==True:
        obj.ham_dia.set(0,N-1, -V)
        obj.ham_dia.set(N-1,0, -V)


    Epot = 0.0
    for i in xrange(N):
        x = q.get(i);  Epot += x*x
    Epot = 0.5*k*Epot

    for i in xrange(N):
        x = q.get(i);
        obj.ham_dia.set(i,i, Epot + alpha*x*(1.0+0.0j))


    for i in xrange(N):
        for n in xrange(N):
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
   
