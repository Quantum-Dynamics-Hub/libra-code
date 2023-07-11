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
###################################################################
#
# This example demonstrates how to obtain normal modes from MD data
# 
###################################################################

import os
import math
import sys
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

def run_test(syst):

    params = {"nsteps":100, "nperiods":1, "scale":250.0}

    if syst=="Si8":
        print "================ Si8 =============="         

        R, V, A, M, E = QE_methods.read_md_data("x0.xml")

        params.update({"print_modes":[5]})

        params.update({"prefix":"_Si8_cov_flag0", "cov_flag":0 })
        res = normal_modes.compute_cov( R, V, A, M, E, params)
        # res is a tuple (w, w_inv_cm, U_v,  w2, w2_inv_cm, U_a)

        params.update({"prefix":"_Si8_cov_flag1", "cov_flag":1 })
        res = normal_modes.compute_cov( R, V, A, M, E, params)
        # res is a tuple (w, w_inv_cm, U_v,  w2, w2_inv_cm, U_a)

        T = 300.0
        params.update({"prefix":"_Si8_cov2_flag0", "cov_flag":0 })
        res = normal_modes.compute_cov2( R, A, M, E, T, params)
        # res is a tuple (w_a, w_inv_cm, U_a)

        params.update({"prefix":"_Si8_cov2_flag1", "cov_flag":1 })
        res = normal_modes.compute_cov2( R, A, M, E, T, params)
        # res is a tuple (w_a, w_inv_cm, U_a)

    elif syst=="H2":
        print "================ H2 =============="
        R, V, A, M, E = QE_methods.read_md_data("x0_h2.xml")

        params.update({"print_modes":[4]})

        params.update({"prefix":"_H2_cov_flag0", "cov_flag":0 })
        res = normal_modes.compute_cov( R, V, A, M, E, params)
        # res is a tuple (w, w_inv_cm, U_v,  w2, w2_inv_cm, U_a)

        params.update({"prefix":"_H2_cov_flag1", "cov_flag":1 })
        res = normal_modes.compute_cov( R, V, A, M, E, params)                                       
        # res is a tuple (w, w_inv_cm, U_v,  w2, w2_inv_cm, U_a)

        T = 300.0        
        params.update({"prefix":"_H2_cov2_flag0", "cov_flag":0 })
        res = normal_modes.compute_cov2( R, A, M, E, T, params)
        # res is a tuple (w_a, w_inv_cm, U_a)

        params.update({"prefix":"_H2_cov2_flag1", "cov_flag":1 })
        res = normal_modes.compute_cov2( R, A, M, E, T, params)
        # res is a tuple (w_a, w_inv_cm, U_a)


    elif syst=="CO2":
        print "================ CO2 =============="
        R, V, A, M, E = QE_methods.read_md_data("x0_co2.xml")

        params.update({"print_modes":[4]})

        params.update({"prefix":"_CO2_cov_flag0", "cov_flag":0 })
        res = normal_modes.compute_cov( R, V, A, M, E, params)
        # res is a tuple (w, w_inv_cm, U_v,  w2, w2_inv_cm, U_a)

        params.update({"prefix":"_CO2_cov_flag1", "cov_flag":1 })
        res = normal_modes.compute_cov( R, V, A, M, E, params)                                       
        # res is a tuple (w, w_inv_cm, U_v,  w2, w2_inv_cm, U_a)

        T = 300.0
        params.update({"prefix":"_CO2_cov2_flag0", "cov_flag":0 })
        res = normal_modes.compute_cov2( R, A, M, E, T, params)
        # res is a tuple (w_a, w_inv_cm, U_a)

        params.update({"prefix":"_CO2_cov2_flag1", "cov_flag":1 })
        res = normal_modes.compute_cov2( R, A, M, E, T, params)
        # res is a tuple (w_a, w_inv_cm, U_a)


run_test("Si8")
run_test("H2")
run_test("CO2")
