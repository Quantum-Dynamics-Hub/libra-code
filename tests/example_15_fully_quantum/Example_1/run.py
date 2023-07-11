#*********************************************************************************
#* Copyright (C) 2018  Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
 This tutorial illustrates how to prepare a Gaussian wavefunction and compute its 
 properties

"""


import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *



wfc = Wfcgrid(-15.0, 25.0, 0.01, 1)  # Define a grid [-15, 25] with the dx = 0.01 and 1 electronic state
wfc.init_wfc_1D(-0.2, 5.0, 1.0, 0)   # Initialize a wavefunction at x = -0.2 and dx = 1.0 with initial momentum px = 5.0 
                                     # it starts on the lowevest-level state (0). Well, for this example this is the 
                                     # only state

# Print out some properties of the wavefunction
print "norm = ", wfc.norm()
print "x0 = ", wfc.get_x_1D()
print "px0 = ", wfc.get_px_1D()
print "E_kin = ", wfc.e_kin_1D(2000.0)


# Print out the function in real and reciprocal spaces manually
f = open("_out.txt", "w")
f.close()

for nx in xrange(wfc.Nx):
    x = wfc.xmin + nx*wfc.dx
    kx = wfc.kxmin + nx/(wfc.Nx*wfc.dx)   # reciprocal space vectors ("wavevectors") - are used internally = k-points grid
    px = 2.0*math.pi*kx                   # physical space momentum vectors - directly related to kinetic energy

    f = open("_out.txt", "a")
    f.write("%8.5f %8.5f  %8.5f %8.5f  %8.5f %8.5f\n" % (x, px, wfc.PSI[0].get(nx,0).real, wfc.PSI[0].get(nx,0).imag, 
            wfc.reciPSI[0].get(nx,0).real,  wfc.reciPSI[0].get(nx,0).imag) )
    f.close()


# Print out the function in real and reciprocal spaces manually
i = 0
wfc.print_wfc_1D("_wfc", i, 0)             # print out real-space wavefunction
wfc.print_reci_wfc_1D("_reci_wfc", i, 0)   # print out reciprocal-space wavefunction

