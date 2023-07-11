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


#plot_pes(params)