#*********************************************************************************
#* Copyright (C) 2018 Alexey V. Akimov
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


# Fisrt, we add the location of the library to test to the PYTHON path
#if sys.platform=="cygwin":
#    from cyglibra_core import *
#elif sys.platform=="linux" or sys.platform=="linux2":
#    from liblibra_core import *
#from libra_py import *


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

params = doubleList()

# Electronic - 2 levels, starting at 0-th state
istate = 0
Coeff = CMATRIX(2,1)
Coeff.set(0, 1.0+0.0j)
Coeff.set(1, 0.0+0.0j)

# Nuclear DOFs
mol = Nuclear(1)
mol.mass[0] = 2000.0
mol.q[0] = -10.0
mol.p[0] = 20.0


f = open("_tsh.txt","w")
dt = 1.0

# Initialization
model_SAC(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
ham.compute_adiabatic(1)

epot = Hadi.get(istate, istate).real
mol.f[0] = -d1ham_adi[0].get(istate, istate).real
Hvib = Hadi - 1.0j * dc1_adi[0] * mol.p[0]/mol.mass[0]


use_boltz_factor = 0
T = 300.0
rnd = Random()

for i in xrange(2500):

    # el-dyn    
    propagate_electronic(0.5*dt, Coeff, Hvib)

    # Nuclear dynamics
    mol.propagate_p(0.5*dt)
    mol.propagate_q(dt)

    model_SAC(Hdia, Sdia, d1ham_dia, dc1_dia, mol.q, params)
    ham.compute_adiabatic(1)

    epot = Hadi.get(istate, istate).real
    mol.f[0] = -d1ham_adi[0].get(istate, istate).real
   

    mol.propagate_p(0.5*dt)

    # el-dyn
    Hvib = Hadi - 1.0j * dc1_adi[0] * mol.p[0]/mol.mass[0]
    propagate_electronic(0.5*dt, Coeff, Hvib)


    ekin = compute_kinetic_energy(mol)


    # Now, incorporate surface hop

    # Just choose the TSH scheme below
    g = compute_hopping_probabilities_fssh(Coeff, Hvib, dt)

 
    do_reverse = 1
    ksi = rnd.uniform(0.0, 1.0)

    #print i, ksi
    p_old = mol.p[0]
    f_old = mol.f[0]

    istate_old = istate
    istate_new = hop(istate_old, g, ksi)
    istate = istate_new
   
#    istate = rescale_velocities_adiabatic(mol.p, mol.mass, Hadi, dc1_adi, istate_new, istate_old, do_reverse)

    ekin = compute_kinetic_energy(mol)
    epot = Hadi.get(istate, istate).real
    mol.f[0] = -d1ham_adi[0].get(istate, istate).real    


#    if istate is not istate_old:
#        print "old p = ", p_old, "Ekin_old = ", 0.5 *p_old**2 / mol.mass[0]  , "Epot_old = ", ham.H(state_old,state_old).real, "etot_old = ", 0.5 *p_old**2 / mol.mass[0] + ham.H(state_old,state_old).real
#        print "new p = ",mol.p[0], "Ekin_new = ", 0.5 *mol.p[0]**2 / mol.mass[0], "Epot_new = ", ham.H(el.istate,el.istate).real, "etot_new = ", 0.5 *mol.p[0]**2 / mol.mass[0] + ham.H(el.istate,el.istate).real
#        print "ekin = ", ekin, " epot = ", epot
#        print "f_old = ", f_old, "f_new = ", mol.f[0], " dHdR = ", ham.dHdq(el.istate,el.istate, 0)


    print i, ekin, epot, ekin+epot, istate, mol.f[0], -d1ham_adi[0].get(0,0).real, -d1ham_adi[0].get(1,1).real

    Denmat = Coeff * Coeff.H()

    f.write("i= %3i q[0]= %8.5f p[0]= %8.5f  ekin= %8.5f  epot= %8.5f  etot= %8.5f  |c0|^2= %8.5f  |c1|^2= %8.5f  Re|c01|= %8.5f istate= %8.5f\n" % (i, mol.q[0], mol.p[0], ekin, epot, ekin+epot, Denmat.get(0,0).real, Denmat.get(1,1).real, Denmat.get(0,1).real, istate))

      
f.close()

