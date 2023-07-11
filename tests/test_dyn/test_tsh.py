#*********************************************************************************
#* Copyright (C) 2018 Brendan A. Smith, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
   This file replaces the test_tsh.py  and we use the nHamiltonian object

"""

import os
import sys
import math

"""
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/math_random")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cygrandom import *
    from cygdyn import *
    from cyghamiltonian import *
    from cyglinalg import *
"""

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


# First, create the Hamiltonian
ham = nHamiltonian(2,2,1)  

# Allocate memory and connect the allocated objects to the ham object
Hdia = CMATRIX(2,2);   Sdia = CMATRIX(2,2);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
invSdia = CMATRIX(2,2);
Hadi = CMATRIX(2,2);   ham.set_ham_adi_by_ref(Hadi);  
U = CMATRIX(2,2);      ham.set_basis_transform_by_ref(U); 
Cdia = CMATRIX(2,1);   ham.set_ampl_dia_by_ref(Cdia)
Cadi = CMATRIX(2,1);   ham.set_ampl_adi_by_ref(Cadi)

d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
for i in [0]:
    d1ham_dia.append( CMATRIX(2,2) ); d1ham_adi.append( CMATRIX(2,2) )
    dc1_dia.append( CMATRIX(2,2) );   dc1_adi.append( CMATRIX(2,2) )

ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)

# Electronic - 2 levels, starting at 1-st state
istate = 1
Coeff = CMATRIX(2,1)
Coeff.set(0, 0.0+0.0j)
Coeff.set(1, 1.0+0.0j)

# Nuclear DOFs
mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = -10.0
mol.p[0] = 20.0

f = open("_tsh.txt","w")
dt = 10.0


# Initialization
params = set_params_2S_1D_sin("SIN:Case3:JCC:2016:37:1626")
print Cpp2Py(params)

model_2S_1D_sin(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
ham.compute_adiabatic(1)

epot = Hadi.get(istate, istate).real
mol.f[0] = -d1ham_adi[0].get(istate, istate).real
Hvib = Hadi - 1.0j * dc1_adi[0] * mol.p[0]/mol.mass[0]


use_boltz_factor = 0
T = 300.0
rnd = Random()


for i in xrange(5000):

    # el-dyn    
    propagate_electronic(0.5*dt, Coeff, Hvib)

    # Nuclear dynamics
    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)

    model_2S_1D_sin(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
    ham.compute_adiabatic(1)

    epot = Hadi.get(istate, istate).real
    mol.f[0] = -d1ham_adi[0].get(istate, istate).real

    mol.propagate_p(0.5*dt)
    # End Nuclear dynamics

    # el-dyn
    Hvib = Hadi - 1.0j * dc1_adi[0] * mol.p[0]/mol.mass[0]
    propagate_electronic(0.5*dt, Coeff, Hvib)

    ekin = compute_kinetic_energy(mol)



    # Now, incorporate surface hop
    # Just choose the TSH scheme below
    g = compute_hopping_probabilities_fssh(Coeff, Hvib, dt)

    do_reverse = 1
    ksi = rnd.uniform(0.0, 1.0)

    istate_old = istate                   # the current state
    istate_new = hop(istate_old, g, ksi)  # proposed state
    istate = rescale_velocities_adiabatic(mol.p, mol.mass, Hadi, dc1_adi, istate_new, istate_old, do_reverse)


    ekin = compute_kinetic_energy(mol)
    epot = Hadi.get(istate, istate).real
    mol.f[0] = -d1ham_adi[0].get(istate, istate).real    
    Denmat = Coeff * Coeff.H()

    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f  |c0|^2= %8.5f  |c1|^2= %8.5f  Re|c01|= %8.5f istate= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot, Denmat.get(0,0).real, Denmat.get(1,1).real, Denmat.get(0,1).real, istate))

      
f.close()

