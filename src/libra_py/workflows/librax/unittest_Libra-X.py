#*********************************************************************************
#* Copyright (C) 2017 Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file unittest_libra-X.py
# This script is for testing different functions in Libra-X
# using PYTHON unittest module

import os
import sys
import math


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from hamiltonian_vib import *
import unittest


########################################################
#                    UNIT TEST
########################################################

class Unittest_Force_orthogonal(unittest.TestCase):
    def test_force_orthogonal(self):
    ##
    # This function tests if firce_orthogonal giving desired output.
    # Here we input arbitrary overlap matrix, orthogonal transformation matrix, and
    # atomic gradients and resulting output is tested against desired value.

        # First generate sample inputs for unittest: smat, cmt and all_grads
        # A simple 2 state 2 atomic ststem
        # Consider S_01 = 0.5, then the smat becomes,
        smat=CMATRIX(2,2)
        smat.set(0,0,1.0)
        smat.set(1,1,1.0)
        smat.set(0,1,0.5)
        smat.set(1,0,0.5)

        # Construct an orthogonal transformation matrix as the following
        cmt =CMATRIX(2,2)
        cmt.set(0,0,0.8)
        cmt.set(1,1,0.8)
        cmt.set(0,1,0.6)
        cmt.set(1,0,0.6)

        # Now create atomic froce components
        all_grads=[]
        for i in xrange(2):  # Number of electronic states
            grads=[]
            for j in xrange(2):  # number of atoms
                g=VECTOR(0.1+i*0.01,-0.2+i*0.01,0.3+i*0.01)
                grads.append(g)
            all_grads.append(grads)

        # The expected orthogonal force F[0][i].(x,y,z) = (22.0*f[0][i].(x,y,z)+15.0*f[1][i].(x,y,z)/25.0
        #                           and F[1][i].(x,y,z) = (15.0*f[0][i].(x,y,z)+22.0*f[1][i].(x,y,z)/25.0, 
        # where f being non-orthogonal forces
        # Example tests are given below
        self.assertAlmostEqual( force_orthogonal(smat,cmt,all_grads)[0][0].x,0.04*(22.0*all_grads[0][0].x+15.0*all_grads[1][0].x))
        #self.assertAlmostEqual( force_orthogonal(smat,cmt,all_grads)[0][0].y,0.04*(22.0*all_grads[0][0].y+15.0*all_grads[1][0].y))
        #self.assertAlmostEqual( force_orthogonal(smat,cmt,all_grads)[0][0].z,0.04*(22.0*all_grads[0][0].z+15.0*all_grads[1][0].z))
        # This tests accuracy up to 7 decimal places

if __name__ =="__main__":
    #smat,cmt,all_grads=input_force_orthogonal()
    unittest.main()

