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

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



print "\nTest 2: Kinetic and nuclear energy of normalized Gaussians"
f = open("3D_eri.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(1.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

for i in range(0,50):
    Rb.x = 0.1 * i

    sasa_sbsb   = electron_repulsion_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  0,0,0, 1.3, Rb)
    sasa_sasb   = electron_repulsion_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Ra, 0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb)
    sasb_sasb   = electron_repulsion_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, 0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb)

    sapxa_sbpxb = electron_repulsion_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  1,0,0, 1.3, Rb)
    sapxa_sapxb = electron_repulsion_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    sapxa_sbpxa = electron_repulsion_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  1,0,0, 1.3, Ra)
    pxapxa_sbpxa = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 0,0,0, 1.3, Rb,  1,0,0, 1.3, Ra)


    pxapxa_pxbpxb   = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 1,0,0, 1.3, Rb,  1,0,0, 1.3, Rb)
    pxapxa_pxapxb   = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Ra, 1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    pxapxb_pxapxb   = electron_repulsion_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, 1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (Rb.x, sasa_sbsb, sasa_sasb, sasb_sasb, sapxa_sbpxb, sapxa_sapxb, sapxa_sbpxa, pxapxa_sbpxa, pxapxa_pxbpxb, pxapxa_pxapxb, pxapxb_pxapxb) )


f.close()



