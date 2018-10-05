#*********************************************************************************                     
#* Copyright (C) 2018 Alexey V. Akimov                                                   
#*                                                                                                     
#* This file is distributed under the terms of the GNU General Public License                          
#* as published by the Free Software Foundation, either version 2 of                                   
#* the License, or (at your option) any later version.                                                 
#* See the file LICENSE in the root directory of this distribution   
#* or <http://www.gnu.org/licenses/>.          
#***********************************************************************************
## \file models_Holstein.py 
#
# This module implements Holstein Hamiltonians
#
#
import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

class tmp:
    pass    


def Holstein_uncoupled(q, params):
    """
    Implementation of a generic Holstein Hamiltonian

    \param[in] q [ndof x 1, MATRIX] coordinate of the particle
    \param[in] params [dictionary] parameters of the model
 
    """

    ndof = q.num_of_rows  # the number of nuclear DOFs
    N = ndof              # in this case, each site has one nuclear DOF

    k = params["k_harmonic"]  # force constant
    alpha = params["el-phon_coupling"] # local electron-phonon coupling
    V = params["site_coupling"]  # diabatic coupling between the adjacent sites


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

    if params["is_periodic"]==1:
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
                obj[n].ham_dia.set(i,i, (k*x + alpha)*(1.0+0.0j))  #  dH(i,i)/dx_n
            else:
                obj[n].ham_dia.set(i,i, k*x*(1.0+0.0j))            #  dH(i,i)/dx_n

    return obj



def get_Holstein_set1():
    """
    Parameters from:
    Qiu, J.; Bai, X.; Wang, L. Crossing Classified and Corrected Fewest Switches Surface Hopping.
    J. Phys. Chem. Lett. 2018, 9, 4319-4325.
    """

    cm1 = 4.5563e-6  # cm^-1 in Ha units
    Ang = 1.0/0.529177249  #  Angstrom in Bohr units
    amu = 1822.89  # amu in units of electron mass
    ps = 4.134137e4  # ps in atomic time units

    params = {}
    params["k_harmonic"]  = 14500 * amu/(ps*ps),
    params["el-phon_coupling"] = 3500.0 * (cm1/Ang) # in Ha/Bohr
    params["mass"] = 250.0 * amu
    params["site_coupling"] = 10.0 * cm1
    params["is_periodic"] = 0

    return params
   
