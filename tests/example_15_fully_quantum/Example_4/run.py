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

import harmonic

class tmp:
    pass


def model1(q, params, full_id):
    """
    Harmonic potential
    Hdia = 0.5*k*(x-x0)^2   
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
    x0, k = params["x0"], params["k"]

    Hdia.set(0,0, 0.5*k*(x-x0)*(x-x0)*(1.0+0.0j) )
    Sdia.set(0,0, 1.0+0.0j)

    for i in [0]:
        #  d Hdia / dR_0
        d1ham_dia[i].set(0,0, k*(x-x0)*(1.0+0.0j) )

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

    x0 = -50
    q = MATRIX(1,1)

    f = open("_pes.txt", "w")
    f.close()

    for i in xrange(1000):
        x = x0 + 0.1*i
        q.set(0,0,x)
        full_id = Py2Cpp_int([0,0])
        obj = compute_model(q, params, full_id)

        f = open("_pes.txt", "a")
        f.write( "%8.5f   %8.5f   %8.5f \n" % (x, obj.ham_dia.get(0,0).real, obj.d1ham_dia[0].get(0,0).real))
        f.close()






def run_exact(params, case):
    """
    The main routine to run fully quantum calculations
    """

    # Here we initialize the grid and wavefunction
    wfc = Wfcgrid(-20.0, 20.0, 0.01, 1)

    wfc.init_wfc_1D_HO(Py2Cpp_int(params["wfc"]["init_state"]), 
                       Py2Cpp_int(params["wfc"]["nu"]),
                       Py2Cpp_complex(params["wfc"]["weights"]),
                       Py2Cpp_double(params["wfc"]["x0"]),
                       Py2Cpp_double(params["wfc"]["px0"]),
                       Py2Cpp_double(params["wfc"]["alpha"]) ) 

#    wfc.normalize_wfc_1D()

    wfc.update_potential_1D(potential, params)
    dt = params["dt"]
    wfc.update_propagator_1D(0.5*dt, params["mass"])  # this is important because we are using exp(-0.5*dt*H_loc)...
    wfc.update_propagator_K_1D(dt, params["mass"])    # ... together with exp(-dt*H_non-loc)            

    f = open("_pops"+str(case)+".txt", "w")
    f.close()
    os.system("mkdir _res"+str(case)+"")

    # Compute dynamics
    cum = 0.0

    exp_pow = doubleList()
    exp_pow.append(2.0)

    for i in xrange(params["nsnaps"]):  # time steps

        epot = wfc.e_pot_1D()
        ekin = wfc.e_kin_1D(params["mass"])
        etot = epot + ekin
        x,  px  = wfc.get_x_1D(), wfc.get_px_1D()
        x2  = wfc.get_pow_x_1D(exp_pow[0])
        px2 = wfc.get_pow_px_1D(2)

        f = open("_pops"+str(case)+".txt", "a")
        f.write("%8.5f   %8.5f    %8.5f   %8.5f  %8.5f  %8.5f    %8.5f    %8.5f  %8.5f  %8.5f\n"
                % (i*params["nsteps"]*dt, ekin, epot, etot, cum, x, px, x2, px2, wfc.norm_1D()))
        f.close()

        wfc.print_wfc_1D("_res"+str(case)+"/wfc", i, 0)   # Be sure to make a "res" directory


        for j in xrange(params["nsteps"]):  # time steps
            wfc.propagate_exact_1D(0)

            res = Py2Cpp_double([0.0])
            wfc.flux_1D(params["barrier"], res, params["mass"])
            cum = cum + res[0]*dt

        

def test(case):
    """
    This function intiializes the params dictionary for a given case    
    """

    params = { "mass":2000.0, "dt":1.0, "barrier":0.00, "nsnaps":200, "nsteps":10 }
    params.update( {"x0":0.0, "k":0.032} )

    # Compute omega based on k
    omega = math.sqrt(params["k"]/params["mass"])
    params.update({"omega":omega})

    if case == 0:

        wfc = {}
        wfc.update({"init_state":[0]})
        wfc.update({"nu":[0]})
        wfc.update({"weights":[1.0+0.0j]})
        wfc.update({"x0":[0.0]})
        wfc.update({"px0":[0.0]})
        wfc.update({"alpha":[ math.sqrt(params["k"] * params["mass"]) ]})

        params.update( {"model":1, "wfc": wfc} )

    elif case == 1:

        wfc = {}
        wfc.update({"init_state":[0, 0]})
        wfc.update({"nu":[0, 1]})
        wfc.update({"weights":[1.0+0.0j, 1.0+0.0j]})
        wfc.update({"x0":[0.0, 0.0]})
        wfc.update({"px0":[0.0, 0.0]})
        alp = math.sqrt(params["k"] * params["mass"])
        wfc.update({"alpha":[ alp, alp ] })

        params.update( {"model":1, "wfc": wfc} )

    elif case == 2:

        wfc = {}
        wfc.update({"init_state":[0, 0, 0, 0]})
        wfc.update({"nu":[0, 1, 2, 3]})
        wfc.update({"weights":[1.0+0.0j, -1.0+0.0j, -1.0+0.0j, 1.0+0.0j]})
        wfc.update({"x0":[0.0, 0.0, 0.0, 0.0]})
        wfc.update({"px0":[0.0, 0.0, 0.0, 0.0]})
        alp = math.sqrt(params["k"] * params["mass"])
        wfc.update({"alpha":[ alp, alp, alp, alp ] })

        params.update( {"model":1, "wfc": wfc} )


    harmonic.run_analytical(params)
    run_exact(params, case)

    plot_pes(params)


test(1)


# k = m*w^2 => w = sqrt(k/m)
# alpha = m * w = m * sqrt(k/m) = sqrt(m*k)
