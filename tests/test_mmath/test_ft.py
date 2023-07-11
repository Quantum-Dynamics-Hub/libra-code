#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
   This files demonstrates the 1D and 2D discrete (FT) and fast (FFT) Fourier transform
   both direct and inverse   

"""

import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def test_1D():

    rnd = Random()
    N = int(math.pow(2, 10))

    print "N = ", N

    x = CMATRIX(N, 1); y = CMATRIX(N, 1); z = CMATRIX(N, 1)

    for i in xrange(N):
        x.set(i, 0, rnd.normal(), rnd.normal() )



    print "============= Doing the DFT and inverse DFT =================="
    t = Timer(); t.start()

    dft(x, y)
    inv_dft(y, z)

    t.stop(); print "Time spent = ", t.show(), " s"
    err = z - x; print (err.H() * err).get(0)


    print "============= Doing the 1D-CFT (version 1) and inverse 1D-CFT (version 1) =================="
    t = Timer(); t.start()

    xmin, dx = -10.0, 0.1
    cft(x, y, xmin, dx)
    inv_cft(y, z, xmin, dx)

    t.stop(); print "Time spent = ", t.show(), " s"
    err = z - x; print (err.H() * err).get(0)


    print "============= Doing the 1D-CFT (version 2) and inverse 1D-CFT (version 2) =================="
    t = Timer(); t.start()

    xmin, kmin, dx = -10.0, -10.0, 0.1
    cft(x, y, xmin, kmin, dx)
    inv_cft(y, z, xmin, kmin, dx)

    t.stop(); print "Time spent = ", t.show(), " s"
    err = z - x; print (err.H() * err).get(0)


    print "============= Doing the 1D-CFFT (version 1) and inverse 1D-CFFT (version 1) =================="
    t = Timer(); t.start()

    xmin, kmin, dx = -10.0, -10.0, 0.1
    cfft(x, y, xmin, kmin, dx)
    inv_cfft(y, z, xmin, kmin, dx)

    t.stop(); print "Time spent = ", t.show(), " s"
    err = z - x; print (err.H() * err).get(0)



def test_2D():

    rnd = Random()
    N = int(math.pow(2, 13))

    print "N = ", N

    x = CMATRIX(N, 2); y = CMATRIX(N, 2); z = CMATRIX(N, 2)

    for i in xrange(N):
        x.set(i, 0, rnd.normal(), rnd.normal() )
        x.set(i, 1, rnd.normal(), rnd.normal() )



    print "============= Doing the 2D-CFT (version 1) and inverse 2D-CFT (version 1) =================="
    t = Timer(); t.start()

    xmin1, kmin1, dx1 = -10.0, -10.0, 0.1
    xmin2, kmin2, dx2 = -10.0, -10.0, 0.1
    cft_2D(x, y, xmin1, xmin2, kmin1, kmin2, dx1, dx2)
    inv_cft_2D(y, z, xmin1, xmin2, kmin1, kmin2, dx1, dx2)

    t.stop(); print "Time spent = ", t.show(), " s"
    err = z - x; print (err.H() * err).tr()


    print "============= Doing the 2D-CFFT (version 1) and inverse 2D-CFFT (version 1) =================="
    t = Timer(); t.start()

    xmin1, kmin1, dx1 = -10.0, -10.0, 0.1
    xmin2, kmin2, dx2 = -10.0, -10.0, 0.1
    cfft_2D(x, y, xmin1, xmin2, kmin1, kmin2, dx1, dx2)
    inv_cfft_2D(y, z, xmin1, xmin2, kmin1, kmin2, dx1, dx2)

    t.stop(); print "Time spent = ", t.show(), " s"
    err = z - x; print (err.H() * err).tr()



test_1D()
test_2D()

