#*********************************************************************************
#* Copyright (C) 2017-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
  This file demonstrates how to run an single Ehrenfest trajectory calculation
  using the built-in Ehrenfest0 function
 
"""

import cmath
import math
import os
import sys


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


"""
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/nHamiltonian_Generic")
sys.path.insert(1,cwd+"/../../_build/src/dyn")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygdyn import *
    from cygconverters import *
    from cygnhamiltonian_generic import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libdyn import *
    from libconverters import *
    from libnhamiltonian_generic import *
    from liblinalg import *
"""

from libra_py import *



def compute_model(q, params, full_id):

    model = params["model"]
    res = None

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    if model==1:
        res = models_Libra.model1(q.col(indx), params)
    elif model==2:
        res = models_Libra.model2(q.col(indx), params)
    elif model==3:
        res = models_Libra.model3(q.col(indx), params)
    elif model==4:
        res = models_Libra.model4(q.col(indx), params)

    res.rep = params["rep"]    
    return res
    


def compute_etot(ham, p, Cdia, Cadi, iM, rep):
    Etot = 0.0
    if rep==0:
        Etot = 0.5*(p.T()*iM*p).get(0) + ham.Ehrenfest_energy_dia(Cdia).real 
    
    elif rep==1:
        Etot = 0.5*(p.T()*iM*p).get(0) + ham.Ehrenfest_energy_adi(Cadi).real        

    return Etot





def run_test(model, rep, outname):
    """
    model - setup the Hamiltonian
    rep - 0 - diabatic, 1 - adiabatic
    outname - the name of the output file
    """

    ham = nHamiltonian(2,2,1)  

    # Allocate memory
    Hdia = CMATRIX(2,2);   Sdia = CMATRIX(2,2);    ham.set_ham_dia_by_ref(Hdia);    ham.set_ovlp_dia_by_ref(Sdia); 
    NACdia = CMATRIX(2,2); NACadi = CMATRIX(2,2);  ham.set_nac_dia_by_ref(NACdia);  ham.set_nac_adi_by_ref(NACadi);
    Hvibdia = CMATRIX(2,2); Hvibadi = CMATRIX(2,2);  ham.set_hvib_dia_by_ref(Hvibdia);  ham.set_hvib_adi_by_ref(Hvibadi);
    invSdia = CMATRIX(2,2);
    Hadi = CMATRIX(2,2);   ham.set_ham_adi_by_ref(Hadi);  
    U = CMATRIX(2,2);      ham.set_basis_transform_by_ref(U); 

    d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
    dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
    for i in [0]:
        d1ham_dia.append( CMATRIX(2,2) ); d1ham_adi.append( CMATRIX(2,2) )
        dc1_dia.append( CMATRIX(2,2) );   dc1_adi.append( CMATRIX(2,2) )
    
    ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
    ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)


    #  Set up the models and compute internal variables
    # Initialization
    # Model parameters 
    params = {}
    params["x0"], params["k"], params["D"], params["V"] = 1.0, 0.1, -0.1, 0.05
    params["omega"] = 0.25
    params["model"] = model
    params["rep"] = rep


    # Simulation parameters
    dt = 1.0

    # Dynamical variables and system-specific properties
    Cdia = CMATRIX(2,1);  
    Cadi = CMATRIX(2,1);  
    Cadi.set(0, 1.0+0.0j); Cadi.set(1, 1.0+0.0j);   Cadi *= (1.0/math.sqrt(2.0))  

    q = MATRIX(1,1); q.set(0, 0.1)
    p = MATRIX(1,1); p.set(0, 0.0)
    iM = MATRIX(1,1); iM.set(0, 1.0/100.0)


    # Initial calculations
    ham.compute_diabatic(compute_model, q, params)
    ham.compute_adiabatic(1); 
    ham.ampl_adi2dia(Cdia, Cadi)

    if rep==0:
        ham.compute_nac_dia(p, iM);  
        ham.compute_hvib_dia(); 
    elif rep==1:
        ham.compute_nac_adi(p, iM);
        ham.compute_hvib_adi(); 

    
    Etot = compute_etot(ham, p, Cdia, Cadi, iM, rep)

    out = open(outname, "w")
    out.close()


    # Do the propagation
    for i in xrange(500):

        if rep==0:
            Ehrenfest0(dt, q, p, iM, Cdia, ham, compute_model, params, 0)
        elif rep==1:
            Ehrenfest0(dt, q, p, iM, Cadi, ham, compute_model, params, 1)

     
        #=========== Properties ==========
        if rep==0:
            ham.ampl_dia2adi(Cdia, Cadi)
        elif rep==1:
            ham.ampl_adi2dia(Cdia, Cadi)



        dm_dia, dm_adi = tsh.compute_dm(ham, Cdia, Cadi, rep, 0)
        Etot = compute_etot(ham, p, Cdia, Cadi, iM, rep)

        out = open(outname, "a")
        ret = (i*dt, q.get(0), p.get(0), Etot, 0.5*(p.T()*iM*p).get(0), dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real,
               Hdia.get(0,0).real, Hdia.get(1,1).real, Hadi.get(0,0).real, Hadi.get(1,1).real, dc1_adi[0].get(0,1).real
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" %  ret )
        out.close()
                                                                    

model = 2

run_test(model, 0, "_0_new.txt")
run_test(model, 1, "_1_new.txt")

        
