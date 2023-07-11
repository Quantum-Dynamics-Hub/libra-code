#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
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
sys.path.insert(1,cwd+"/../_build/src/hamiltonian/Hamiltonian_Atomistic/Hamiltonian_MM")
sys.path.insert(1,cwd+"/../_build/src/converters")
sys.path.insert(1,cwd+"/../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cygconverters import *
    from cyghamiltonian_mm import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from libconverters import *
    from libhamiltonian_mm import *
    from liblinalg import *







class Test_VdW_Interactions(unittest.TestCase):
    """ Summary of the tests:
    """

    def test_1(self):
        """Set coordinates"""
        R = Py2Cpp_VECTOR( [VECTOR(0.0, 0.0, 0.0),  VECTOR(1.0, 0.0, 0.0),  VECTOR(2.0, 0.0, 0.0)  ] )
        i1 = VdW_Interaction()
        i1.set_coords(R, Py2Cpp_int([0,1]))

    def test_2(self):
        """Set translations"""
        T = Py2Cpp_VECTOR( [VECTOR(0.0, 0.0, 0.0),  VECTOR(3.0, 0.0, 0.0) ] )
        i1 = VdW_Interaction()
        i1.set_transl(T, Py2Cpp_int([0,0]))


    def test_3(self):
        """Set forces"""
        F = Py2Cpp_VECTOR( [VECTOR(0.0, 0.0, 0.0),  VECTOR(0.0, 0.0, 0.0), VECTOR(0.0, 0.0, 0.0) ] )
        i1 = VdW_Interaction()
        i1.set_forces(F, Py2Cpp_int([0,1]))

    def test_4(self):
        """Set functional and parameters"""
        i1 = VdW_Interaction()
        i1.set_functional("LJ")
        i1.set_params({"sigma":1.0, "epsilon":0.05, "scale": 1.0, "D": 0.05, "r0":1.0, "alpha":0.5  })


        self.assertAlmostEqual( i1.sigma, 1.0 );
        self.assertAlmostEqual( i1.epsilon, 0.05 );
        self.assertAlmostEqual( i1.scale, 1.0 );
        self.assertAlmostEqual( i1.D, 0.05 );
        self.assertAlmostEqual( i1.r0, 1.0 );
        self.assertAlmostEqual( i1.alpha, 0.5 );

        self.assertAlmostEqual( i1.energy, 0.0 );
        self.assertEqual( i1.int_type, 21 );
        self.assertEqual( i1.functional, 0 );


    def test_5(self):
        """Calculations: Setup only coordinates, no translations, no forces """
        R = Py2Cpp_VECTOR( [VECTOR(0.0, 0.0, 0.0),  VECTOR(1.0, 0.0, 0.0),  VECTOR(2.0, 0.0, 0.0)  ] )
        T = Py2Cpp_VECTOR( [VECTOR(0.0, 0.0, 0.0),  VECTOR(3.0, 0.0, 0.0) ] )
        F = Py2Cpp_VECTOR( [VECTOR(0.0, 0.0, 0.0),  VECTOR(0.0, 0.0, 0.0), VECTOR(0.0, 0.0, 0.0) ] )

        # 0 - 1
        i01 = VdW_Interaction()
        i01.set_functional("LJ")
        i01.set_params({"sigma":1.0, "epsilon":0.05, "scale": 1.0 })
        i01.set_coords(R, Py2Cpp_int([0,1]))

        i01.compute()
        self.assertAlmostEqual( i01.energy, -0.05 );    
        self.assertAlmostEqual( F[0].x, 0.0 );
        self.assertAlmostEqual( F[1].x, 0.0 );

        # 0 - 2 
        i02 = VdW_Interaction()
        i02.set_functional("LJ")
        i02.set_params({"sigma":1.0, "epsilon":0.05, "scale": 1.0 })
        i02.set_coords(R, Py2Cpp_int([0,2]))

        i02.compute()
        self.assertAlmostEqual( i02.energy, 0.05 * ( (0.5)**12 - 2.0 * (0.5)**6 ) );
        self.assertAlmostEqual( F[0].x,  0.0 );   # Note, the forces have to ve zero, because
        self.assertAlmostEqual( F[2].x,  0.0 );   # the storage was not initialized, so the actual 
                                                  # results will not be stored there
       



        

               


if __name__=='__main__':
    unittest.main()




