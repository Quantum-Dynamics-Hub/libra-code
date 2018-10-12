#***********************************************************
# * Copyright (C) 2017-2018 Brendan A. Smith and Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

def compute_hprime_nosoc(es_curr, info, filename):
    """
    This function computes the matrix elements of the dipole operator for the
    case without SOC, and prints them to them to the file specificed by the 
    function parameter "filename"

    \param[in] es_curr A dictionary containing the data for the g-vectors and pw coefficients
    \param[in] info A dictionary containing the basic information regarding the system, such 
               as recip. lattice vectors 
    \param[in] filename This is the name of the output file where the data will be printed
    """  

    coeff = es_curr["Coeff_dia"][0]
    grid  = es_curr["grid"][0]
    b1    = info["b1"] 
    b2    = info["b2"]
    b3    = info["b3"]

    g = MATRIX(len(grid),3)
    for i in xrange(len(grid)):
        g.set(i,0,grid[i].x)
        g.set(i,1,grid[i].y)
        g.set(i,2,grid[i].z)

    reci = MATRIX(3,3)
    reci.set(0,0,b1.x); reci.set(0,1,b1.y); reci.set(0,2,b1.z)
    reci.set(1,0,b2.x); reci.set(1,1,b2.y); reci.set(1,2,b2.z)
    reci.set(2,0,b3.x); reci.set(2,1,b3.y); reci.set(2,2,b3.z)

    Hprime = compute_Hprime(coeff,g,reci)

    Hprime[0].real().show_matrix("%sx_re" % (filename));  Hprime[0].imag().show_matrix("%sx_im" % (filename))    
    Hprime[1].real().show_matrix("%sy_re" % (filename));  Hprime[1].imag().show_matrix("%sy_im" % (filename))    
    Hprime[2].real().show_matrix("%sz_re" % (filename));  Hprime[2].imag().show_matrix("%sz_im" % (filename)) 



