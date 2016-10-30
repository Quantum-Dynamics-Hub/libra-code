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
f = open("3D_kin_nucl.txt","w")
Ra = VECTOR(0.0, 0.0, 0.0)
Rb = VECTOR(1.0, 0.0, 0.0)
Rc = VECTOR(0.0, 0.0, 0.0)

for i in range(0,50):
    Rb.x = 0.1 * i

    s_s   = kinetic_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb)
    s_s   = s_s + nuclear_attraction_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, Ra)
    s_s   = s_s + nuclear_attraction_integral(0,0,0, 1.3, Ra,  0,0,0, 1.3, Rb, Rb)

    s_px   = kinetic_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    s_px   = s_px + nuclear_attraction_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Ra)
    s_px   = s_px + nuclear_attraction_integral(0,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Rb)

    px_px   = kinetic_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb)
    px_px   = px_px + nuclear_attraction_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Ra)
    px_px   = px_px + nuclear_attraction_integral(1,0,0, 1.3, Ra,  1,0,0, 1.3, Rb, Rb)

    px_py   = kinetic_integral(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb)
    px_py   = px_py + nuclear_attraction_integral(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, Ra)
    px_py   = px_py + nuclear_attraction_integral(1,0,0, 1.3, Ra,  0,1,0, 1.3, Rb, Rb)

    f_f   = kinetic_integral(3,0,0, 1.3, Ra,  3,0,0, 1.3, Rb)
    f_f   = f_f + nuclear_attraction_integral(3,0,0, 1.3, Ra,  3,0,0, 1.3, Rb, Ra)
    f_f   = f_f + nuclear_attraction_integral(3,0,0, 1.3, Ra,  3,0,0, 1.3, Rb, Rb)


    f.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \n" % (Rb.x, s_s, s_px, px_px, px_py, f_f ) )


f.close()



