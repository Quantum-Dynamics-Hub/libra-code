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
 This is an extended test suite for validating the CDFT and CFFT and their inverses
"""

import os
import sys
import math
import copy
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



class TestFT(unittest.TestCase):
    def test_1D(self):
        """Tests 1D complex DFT and FFT algorithms"""

        rnd = Random()
        N = int(math.pow(2, 10))
        x = CMATRIX(N, 1); y = CMATRIX(N, 1); z = CMATRIX(N, 1)

        for i in xrange(N):
            x.set(i, 0, rnd.normal(), rnd.normal() )


        #============= Doing the DFT and inverse DFT ==================
        dft(x, y)
        inv_dft(y, z)
        err = z - x         
        self.assertAlmostEqual( abs((err.H() * err).get(0)),  0.0 , 15);


        #============= Doing the 1D-CFT (version 1) and inverse 1D-CFT (version 1) ==================
        xmin, dx = -10.0, 0.1
        cft(x, y, xmin, dx)
        inv_cft(y, z, xmin, dx)
        err = z - x         
        self.assertAlmostEqual( abs((err.H() * err).get(0)),  0.0 , 15);


        #============= Doing the 1D-CFT (version 2) and inverse 1D-CFT (version 2) ==================
        xmin, kmin, dx = -10.0, -10.0, 0.1
        cft(x, y, xmin, kmin, dx)
        inv_cft(y, z, xmin, kmin, dx)
        err = z - x
        self.assertAlmostEqual( abs((err.H() * err).get(0)),  0.0 , 15);


        #============= Doing the 1D-CFFT (version 1) and inverse 1D-CFFT (version 1) ==================
        xmin, kmin, dx = -10.0, -10.0, 0.1
        cfft(x, y, xmin, kmin, dx)
        inv_cfft(y, z, xmin, kmin, dx)
        err = z - x
        self.assertAlmostEqual( abs((err.H() * err).get(0)),  0.0 , 15);


    def test_2D(self):
        """Tests 2D complex DFT and FFT algorithms"""

        rnd = Random()
        N = int(math.pow(2, 10))
        x = CMATRIX(N, 2); y = CMATRIX(N, 2); z = CMATRIX(N, 2)

        for i in xrange(N):
            x.set(i, 0, rnd.normal(), rnd.normal() )
            x.set(i, 1, rnd.normal(), rnd.normal() )

        #============= Doing the 2D-CFT (version 1) and inverse 2D-CFT (version 1) ==================
        xmin1, kmin1, dx1 = -10.0, -10.0, 0.1
        xmin2, kmin2, dx2 = -10.0, -10.0, 0.1
        cft_2D(x, y, xmin1, xmin2, kmin1, kmin2, dx1, dx2)
        inv_cft_2D(y, z, xmin1, xmin2, kmin1, kmin2, dx1, dx2)
        err = z - x
        self.assertAlmostEqual( abs((err.H() * err).get(0)),  0.0 , 15);


        #============= Doing the 2D-CFFT (version 1) and inverse 2D-CFFT (version 1) ==================
        xmin1, kmin1, dx1 = -10.0, -10.0, 0.1
        xmin2, kmin2, dx2 = -10.0, -10.0, 0.1
        cfft_2D(x, y, xmin1, xmin2, kmin1, kmin2, dx1, dx2)
        inv_cfft_2D(y, z, xmin1, xmin2, kmin1, kmin2, dx1, dx2)
        err = z - x
        self.assertAlmostEqual( abs((err.H() * err).get(0)),  0.0 , 15);



if __name__=='__main__':
    unittest.main()



