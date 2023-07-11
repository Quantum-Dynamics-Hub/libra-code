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
import cmath
import math
import os
import sys
import unittest


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



class tmp:
    pass    

def model1(q, params):
    """

    Hdia, Sdia, d1ham_dia, dc1_dia, 

              k*x^2         V
    Hdia =       V        k*(x-x0)^2 + D

    Sdia =  I

    Ddia  = 0.0

    """

    Hdia = CMATRIX(2,2)
    Sdia = CMATRIX(2,2)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(2,2))
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(2,2))


    x = q.get(0)
    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    Sdia.set(0,0, 1.0+0.0j);  Sdia.set(0,1, 0.0+0.0j);
    Sdia.set(1,0, 0.0+0.0j);  Sdia.set(1,1, 1.0+0.0j);

    Hdia.set(0,0, k*x*x*(1.0+0.0j) );   Hdia.set(0,1, V*(1.0+0.0j));
    Hdia.set(1,0, V*(1.0+0.0j));        Hdia.set(1,1, (k*(x-x0)**2 + D)*(1.0+0.0j));

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) );   d1ham_dia[i].set(0,1, 0.0+0.0j);
        d1ham_dia[i].set(1,0, 0.0+0.0j);   d1ham_dia[i].set(1,1,2.0*k*(x-x0)*(1.0+0.0j));

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j);   dc1_dia[i].set(0,1, 0.0+0.0j);
        dc1_dia[i].set(1,0, 0.0+0.0j);   dc1_dia[i].set(1,1, 0.0+0.0j);

#    Hdia.show_matrix()
#    print Hdia, Sdia

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def test():


    ham = nHamiltonian(2,2,1)  

    # Allocate memory
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


    Cadi.set(0, 1.0+0.0j);
    Cadi.set(1, 0.0+0.0j);


    #  Set up the models and compute internal variables
    # Initialization
    # Model parameters 
    params = {}
    params["x0"], params["k"], params["D"], params["V"] = 1.0, 0.01, -0.1, 0.05
    params["omega"] = 0.25


    q = MATRIX(1,1);  q.set(0, -1.0)
    p = MATRIX(1,1);  p.set(0, 0.0)
 
    dt = 1.0
    iM = MATRIX(1,1);  iM.set(0, 20.0)
      

    f = open("verlet.txt","w")    
    f.close()

    for i in xrange(100):
        Verlet(dt, q, p, iM, ham, model1, params)

        x = q.get(0)
        Edia = ham.get_ham_dia().get(0).real
        Eadi = ham.get_ham_adi().get(0).real
        Etot = ham.get_ham_adi().get(0).real + 0.5*(p.T()*iM*p).get(0)
        f = open("verlet.txt","a")            
        f.write("%5i %8.5f %8.5f %8.5f %8.5f\n" % (i, x, Edia, Eadi, Etot) )
        f.close()       

test()






