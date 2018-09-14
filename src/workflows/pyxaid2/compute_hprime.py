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

    g_sz = len(grid)
    npw  = coeff.num_of_rows

    print "g_sz  = ", g_sz
    print "npw   = ", npw 
    
    scl1 = (1.0+0.0j)
    scl2 = (0.0+1.0j)

    # figure out if we are using a completed wavefunction or not
    if g_sz != npw:
        g_sz = min(g_sz,npw)
        is_compl = 1
    else:
        is_compl = 0


    N = coeff.num_of_cols 
    Hx = CMATRIX(N/2,N/2); Hy = CMATRIX(N/2,N/2); Hz = CMATRIX(N/2,N/2);
    for g in xrange(g_sz):

        hx, hy, hz = 0.0, 0.0, 0.0

        for i in xrange(N/2):  
            for j in xrange(N/2):

                tmp  = coeff.H().get(i,g) * coeff.get(g,j)
             
                gx = scl1*(grid[g].x*b1.x + grid[g].y*b2.x + grid[g].z*b3.x)
                gy = scl1*(grid[g].x*b1.y + grid[g].y*b2.y + grid[g].z*b3.y)
                gz = scl1*(grid[g].x*b1.z + grid[g].y*b2.z + grid[g].z*b3.z)

                if(is_compl==0):
                    hx += tmp.real*gx; hy += tmp.real*gy; hz += tmp.real*gz;
                
                if(is_compl==1): 
                    if(g==0):
                        hx += tmp.real*gx; hy += tmp.real*gy; hz += tmp.real*gz;  # This should give zero!  

                    #Now the Hprime_ matrices are purely imaginary, for the case of gamma-symmetry.
                    else:
                        hx +=  2.0*scl2*tmp.real*gx   
                        hy +=  2.0*scl2*tmp.real*gy             
                        hz +=  2.0*scl2*tmp.real*gz

                Hx.add(i,j,hx)
                Hy.add(i,j,hy)
                Hz.add(i,j,hz)

    Hx.real().show_matrix("%sx_re" % (filename));  Hx.imag().show_matrix("%sx_im" % (filename))    
    Hy.real().show_matrix("%sy_re" % (filename));  Hy.imag().show_matrix("%sy_im" % (filename))    
    Hz.real().show_matrix("%sz_re" % (filename));  Hz.imag().show_matrix("%sz_im" % (filename))  


