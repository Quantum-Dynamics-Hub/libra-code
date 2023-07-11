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

class tmp:
    pass

def model1(q, params, full_id):
    """
    Harmonic potential
    Hdia = 0.5*k*x^2   
    Sdia =  1.0
    Ddia  = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )
  
    x = q.col(indx).get(0)
    x0,k,D,V = params["x0"], params["k"], params["D"], params["V"]

    Hdia.set(0,0, k*x*x*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, 2.0*k*x*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj
    
def model2(q, params, full_id):
    """
    Symmetric Double Well Potential
    Hdia = 0.25*q^4 - 0.5*q^2   
    Sdia = 1.0
    Ddia = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )

    x = q.col(indx).get(0)
    x2 = x*x

    Hdia.set(0,0, (0.25*x2*x2 - 0.5*x2)*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, x*(x2 - 1.0)*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj



def model3(q, params, full_id):
    """
    #Params in:  q = list of trajectory positions
    Cubic potential
    Hdia = A*q^2 - B*q^3   
    Sdia = 1.0
    Ddia = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )

    b, ww = 0.2981, 0.01
    m = params["mass"]
    mw2 = m * ww * ww
         
    x = q.col(indx).get(0)
    x2 = x*x

    en = 0.5*mw2*x2 - (b/3.0)*x2*x
    den = mw2*x - b*x2

    Hdia.set(0,0, en*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, den*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)

    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj


def model4(q, params, full_id):
    """
    #Params in:  q = list of trajectory positions
    Cubic potential with switching 
    Hdia = A*q^2 - B*q^3   
    Sdia = 1.0
    Ddia = 0.0
    """

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    Hdia = CMATRIX(1,1)
    Sdia = CMATRIX(1,1)
    d1ham_dia = CMATRIXList();  d1ham_dia.append( CMATRIX(1,1) )
    dc1_dia = CMATRIXList();  dc1_dia.append( CMATRIX(1,1) )

    b, ww = 0.2981, 0.01
    m = params["mass"]
    mw2 = m * ww * ww
         
    x = q.col(indx).get(0)
    x2 = x*x

    fun = 0.5*mw2*x2 - (b/3.0)*x2*x
    dfun = mw2*x - b*x2

    sw = SWITCH(VECTOR(x,0.0,0.0),VECTOR(0.0, 0.0, 0.0), 0.9, 1.5);

    x = 1.5
    x2 = x * x
    fun_const = 0.5*mw2*x2 - (b/3.0)*x2*x
    dfun_const = mw2*x - b*x2

    en = fun*sw[0] + fun_const*(1.0 - sw[0])
    den = dfun *sw[0] + fun*sw[1].x - fun_const*sw[1].x
   

    Hdia.set(0,0, en*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, den*(1.0+0.0j) )

        #  <dia| d/dR_0| dia >
        dc1_dia[i].set(0,0, 0.0+0.0j)


    obj = tmp()
    obj.ham_dia = Hdia
    obj.ovlp_dia = Sdia
    obj.d1ham_dia = d1ham_dia
    obj.dc1_dia = dc1_dia

    return obj




def compute_model(q, params, full_id):
    """
    Generic calculator of the model Hamiltonians
    """

    model = params["model"]
    res = None

    if model==1:
        res = model1(q, params, full_id)
    if model==2:
        res = model2(q, params, full_id)
    if model==3:
        res = model3(q, params, full_id)
    if model==4:
        res = model4(q, params, full_id)

    return res



def potential(q, params):
    """
    Thin wrapper of the model Hamiltonians that can be used in 
    the fully-quantum calculations
    """

    return compute_model(q, params, Py2Cpp_int([0,0]))
    



def plot_pes(params):
    """
    An auxiliary function to compute the PES profiles
    """

    x0 = -0.6
    q = MATRIX(1,1)

    f = open("_pes.txt", "w")
    f.close()

    for i in xrange(2600):
        x = x0 + 0.001*i
        q.set(0,0,x)
        full_id = Py2Cpp_int([0,0])
        obj = compute_model(q, params, full_id)

        f = open("_pes.txt", "a")
        f.write( "%8.5f   %8.5f   %8.5f \n" % (x, obj.ham_dia.get(0,0).real, obj.d1ham_dia[0].get(0,0).real))
        f.close()




def run_exact(params, opt):
    """
    The main routine to run fully quantum calculations
    """

    # Here we initialize the grid and wavefunction
    wfc = Wfcgrid(-15.0, 25.0, 0.01, 1)
    wfc.init_wfc_1D(params["q0"], params["p0"], params["sq0"], 0)

    wfc.update_potential_1D(potential, params)
    dt = params["dt"]
    wfc.update_propagator_1D(0.5*dt, params["mass"])  # this is important because we are using exp(-0.5*dt*H_loc)...
    wfc.update_propagator_K_1D(dt, params["mass"])    # ... together with exp(-dt*H_non-loc)            


    f = open("_pops.txt", "w")
    f.close()

    # Compute dynamics
    cum = 0.0
    for i in xrange(params["nsnaps"]):  # time steps
        for j in xrange(params["nsteps"] ):  # time steps
            wfc.propagate_exact_1D(0)

            if opt==0:
                """
                This is the case with absorbing potential - the total energy of the remaining
                wavefunction is not conserved! Also, the population dynamics may be delayed 
                depending on where the potential is located
                """
                res = wfc.absorb_1D(5.0)
                cum = cum + res[1][0]

            elif opt==1:
                """
                This is the case with population density fluxes. 
                The total population and energy is conserved. 
                The dynamics is not delayed, because we measure exactly what we want (crossing some point)
                """

                res = Py2Cpp_double([0.0])
                wfc.flux_1D(params["barrier"], res, params["mass"])
                cum = cum + res[0]*dt



        epot = wfc.e_pot_1D()
        ekin = wfc.e_kin_1D(params["mass"])
        etot = epot + ekin


        f = open("_pops.txt", "a")
        f.write("%8.5f   %8.5f  %8.5f  %8.5f    %8.5f\n" % (i*params["nsteps"]*dt, ekin, epot, etot, cum))
        f.close()

        wfc.print_wfc_1D("res/wfc", i, 0)   # Be sure to make a "res" directory

        



# Model parameters 
case = 3
params = { "mass":2000.0, "nsnaps":100, "nsteps":2, "dt":10.0 } 

if case==1:
    # Double well potential
    params.update({"model":2, "q0":-1.1, "p0":0.0, "sq0":0.04, "sp0":0.0, "barrier":0.00 })
elif case==2:
    # Cubic, no switching
    params.update({"model":3, "q0":-0.2, "p0":0.0, "sq0":0.1, "sp0":0.0, "barrier":0.67 })
elif case==3:
    # Cubic, with switching
    params.update({"model":4, "q0":-0.2, "p0":0.0, "sq0":0.1, "sp0":5.25, "barrier":0.67 })   


plot_pes(params)

opt = 1 
run_exact(params, opt)
