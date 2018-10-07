#*********************************************************************************
#* Copyright (C) 2017-2018  Brendan A. Smith, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/


import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

import math
import sys

def compute_energy(n,omega):
    """
    Computes the energy for the nth HO vibraitonal mode    
    with frequency omega
   
    omega - Frequency associated with the harmonic potential
    n - The HO eigenstate index

    Returns: The energy for the nth HO eigenstate
    """

    E = omega * ( n + 0.5 )
    return E


def tot_energy(coeff, omega):
    """
    Computes the total energy of a superposition of HOs
   
    coeff = CMATRIX(nstates,1) coefficients of the superposition
    omega (double) the frequency of the HO

    """

    n = coeff.num_of_rows

    E = 0.0
    for i in xrange(n):
        E = E + ( i + 0.5 ) * (coeff.get(i,0).conjugate() * coeff.get(i,0)).real
    E = E * omega

    return E



def compute_q(params,coeff):
    """
    Computes the expecation value < PSI | q | PSI >. Where PSI 
    is a user defined superposition of HO eigenstates

    params - A dictionary containg simulation parameters
    coeff_curr - a CMATRIX HO eigenstate coefficients of size (nbasis,1)

    Returns: < PSI | q | PSI >
    """


    mass  = params["mass"]
    omega = params["omega"]
    nbasis = coeff.num_of_rows

    q = 0.0
    for n in xrange(nbasis-1):
        c1 = coeff.get(n+1,0).conjugate() * coeff.get(n,0)
        q +=  math.sqrt(n+1) * c1.real

    q *= math.sqrt(2.0/(mass*omega))

    return q


def compute_q2(params,coeff_curr,energy,t):
    """
    Computes the expecation value < PSI | q^2 | PSI >. Where PSI 
    is a user defined superposition of HO eigenstates

    params - A dictionary containg simulation parameters
    coeff_curr - a CMATRIX HO eigenstate coefficients of size (nbasis,1)
    energy = a MATRIX of HO eigenenergies of size (nbasis,1)
    t = current time step

    Returns: < PSI | q^2 | PSI >
    """

    mass  = params["mass"]
    omega = params["omega"]
    nbasis = coeff_curr.num_of_rows

    q = 0.0
    # for ( a_+ a_- ) + ( a_- a_+ )
    for n in xrange(nbasis):
        c1 = coeff_curr.get(n,0).conjugate() * coeff_curr.get(n,0)
        q += (2*n+1) * c1.real

    # for ( a_+^2 ) + ( a_-^2 )
    for n in xrange(nbasis - 2):
        c2 = coeff_curr.get(n+2,0).conjugate() * coeff_curr.get(n,0)
        q += 2.0 * math.sqrt(n+1) * math.sqrt(n+2) * c2.real

    q *= 0.5/(mass*omega)

    return q



def compute_p(params,coeff):
    """
    Computes the expecation value < PSI | p | PSI >. Where PSI 
    is a user-defined superposition of HO eigenstates

    params - A dictionary containg simulation parameters
    coeff_curr - a CMATRIX HO eigenstate coefficients of size (nbasis,1)
    energy = a MATRIX of HO eigenenergies of size (nbasis,1)
    t = current time step

    Returns: < PSI | p | PSI >
    """

    mass  = params["mass"]
    omega = params["omega"]
    nbasis = coeff.num_of_rows

    p = 0.0
    for n in xrange(nbasis-1):
        c1 = coeff.get(n+1,0).conjugate() * coeff.get(n,0)
        p -=  math.sqrt(n+1) * c1.imag
 
    p *= math.sqrt(2.0*mass*omega)

    return p



def compute_p2(params,coeff_curr,energy,t):
    """
    Computes the expecation value < PSI | p^2 | PSI >. Where PSI 
    is a user defined superposition of HO eigenstates

    params - A dictionary containg simulation parameters
    coeff_curr - a CMATRIX HO eigenstate coefficients of size (nbasis,1)
    energy = a MATRIX of HO eigenenergies of size (nbasis,1)
    t = current time step

    Returns: < PSI | p^2 | PSI >
    """

    mass  = params["mass"]
    omega = params["omega"]
    nbasis = coeff_curr.num_of_rows

    p = 0.0

    # for ( a_+ a_- ) + ( a_- a_+ )
    for n in xrange(nbasis):
        c1 = coeff_curr.get(n,0).conjugate() * coeff_curr.get(n,0)
        p += (2*n+1) * c1.real

    # for ( a_+^2 ) + ( a_-^2 )
    for n in xrange(nbasis - 2):
        c2 = coeff_curr.get(n+2,0).conjugate() * coeff_curr.get(n,0)
        p -= 2.0 * math.sqrt(n+1) * math.sqrt(n+2) * c2.real

    p *= 0.5*mass*omega

    return p



def run_analytical(params):    
    """
    This is the driver fucntion for the file harmonic.py 

    nsnaps - number of instances to be printed
    nsteps - number of steps between each nsnap
    params - A dictionary containg simulation parameters
    case   - the index for a particular case
    """

    dt = params["dt"]
    coeff = params["coeff"]
    etot = 0.0


    nbasis = len(params["coeff"])

    # Initialize the coeffcient list
    coeff = CMATRIX(nbasis,1)
    energy = MATRIX(nbasis,1)

    for n in xrange(nbasis):
        coeff.set(n,0, params["coeff"][n])
        energy.set(n,0, compute_energy(n,params["omega"]) )

    # Normalize the wavefunction
    norm = (coeff.H()*coeff).get(0).real
    coeff = coeff / math.sqrt(norm)



    f = open("_res.txt", "w");  
    f.close()


    t = 0
    for i in xrange(params["nsnaps"]): 

        # Compute propeties (expectation values, etc)
        etot = tot_energy(coeff, params["omega"]) 
        q  = compute_q(params,coeff)
        p  = compute_p(params,coeff)
        q2 = compute_q2(params,coeff,energy,t)
        p2 = compute_p2(params,coeff,energy,t)
   
        # Print results
        f = open("_res.txt", "a")
        f.write("%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" % (t, etot, q, p, q2, p2) )
        f.close()

        for j in xrange(params["nsteps"]):
            # For each time step, we have to update all coeffcients 
            for n in xrange(nbasis):
                argg = dt*params["omega"]*(n+0.5)
                U_n = math.cos( argg ) - 1.0j*math.sin( argg )
                coeff.scale(n,0, U_n)

            t += dt
            


