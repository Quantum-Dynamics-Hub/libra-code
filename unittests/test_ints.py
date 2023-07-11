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


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *




class TestInts(unittest.TestCase):
    """ Summary of the tests:
    1 - test overlap integrals
    """

    def test_1(self):
        """Test overlaps"""

        f = open("reference_integrals.txt") 
        A = f.readlines()
        f.close()
       
        for a in A:
            inp = a.split()

            a1 = float(inp[0])
            R1 = VECTOR(float(inp[1]), float(inp[2]), float(inp[3])) 
            n1,m1,k1 = int(float(inp[4])), int(float(inp[5])), int(float(inp[6]))

            a2 = float(inp[7])
            R2 = VECTOR(float(inp[8]), float(inp[9]), float(inp[10]))
            n2,m2,k2 = int(float(inp[11])), int(float(inp[12])), int(float(inp[13]))

            ref_ovlp = float(inp[14])


            g1 = PrimitiveG()
            g2 = PrimitiveG()
            g1.init(n1,m1,k1, a1, R1)
            g2.init(n2,m2,k2, a2, R2)

            ovlp = gaussian_overlap(g1,g2,0)           

               
            self.assertAlmostEqual( ovlp, ref_ovlp )
                
        print "Tested ", len(A), "unnormalized overlap integrals"


    def test_2(self):
        """Test normalized overlaps"""

        f = open("reference_integrals.txt")
        A = f.readlines()
        f.close()

        for a in A:
            inp = a.split()

            a1 = float(inp[0])
            R1 = VECTOR(float(inp[1]), float(inp[2]), float(inp[3]))
            n1,m1,k1 = int(float(inp[4])), int(float(inp[5])), int(float(inp[6]))

            a2 = float(inp[7])
            R2 = VECTOR(float(inp[8]), float(inp[9]), float(inp[10]))
            n2,m2,k2 = int(float(inp[11])), int(float(inp[12])), int(float(inp[13]))

            ref_ovlp = float(inp[15])


            g1 = PrimitiveG()
            g2 = PrimitiveG()
            g1.init(n1,m1,k1, a1, R1)
            g2.init(n2,m2,k2, a2, R2)

            ovlp = gaussian_overlap(g1,g2)


            self.assertAlmostEqual( ovlp, ref_ovlp )

        print "Tested ", len(A), "normalized overlap integrals"


    def test_3(self):
        """Test kinetic integrals"""

        f = open("reference_integrals.txt")
        A = f.readlines()
        f.close()

        for a in A:
            inp = a.split()

            a1 = float(inp[0])
            R1 = VECTOR(float(inp[1]), float(inp[2]), float(inp[3]))
            n1,m1,k1 = int(float(inp[4])), int(float(inp[5])), int(float(inp[6]))

            a2 = float(inp[7])
            R2 = VECTOR(float(inp[8]), float(inp[9]), float(inp[10]))
            n2,m2,k2 = int(float(inp[11])), int(float(inp[12])), int(float(inp[13]))

            ref_kin = float(inp[16])


            g1 = PrimitiveG()
            g2 = PrimitiveG()
            g1.init(n1,m1,k1, a1, R1)
            g2.init(n2,m2,k2, a2, R2)

            kin = kinetic_integral(g1,g2)


            self.assertAlmostEqual( kin, ref_kin )

        print "Tested ", len(A), "kinetic energy integrals"


       



if __name__=='__main__':
    unittest.main()




