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

def compute_hprime(es_curr, info, filename):

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
    Hx = CMATRIX(N,N); Hy = CMATRIX(N,N); Hz = CMATRIX(N,N);
    for i in xrange(N):  
        for j in xrange(N):

            hx, hy, hz = 0.0, 0.0, 0.0
            for g in xrange(g_sz):

                tmp  = coeff.H().row(i).col(g) * coeff.col(j).row(g)

                gx = scl1*(grid[g].x*b1.x + grid[g].y*b2.x + grid[g].z*b3.x)
                gy = scl1*(grid[g].x*b1.y + grid[g].y*b2.y + grid[g].z*b3.y)
                gz = scl1*(grid[g].x*b1.z + grid[g].y*b2.z + grid[g].z*b3.z)

                if(is_compl==0):
                    hx += tmp.get(0)*gx; hy += tmp.get(0)*gy; hz += tmp.get(0)*gz;
                
                if(is_compl==1): 
                    if(g==0):
                        hx += tmp.get(0)*gx; hy += tmp.get(0)*gy; hz += tmp.get(0)*gz;  # This should give zero!  
                    #Now the Hprime_ matrices are purely imaginary, for the case of gamma-symmetry.
                    else:
                        hx +=  2.0*scl2*tmp.get(0).real*gx   
                        hy +=  2.0*scl2*tmp.get(0).real*gy             
                        hz +=  2.0*scl2*tmp.get(0).real*gz

            Hx.set(i,j,hx)
            Hy.set(i,j,hy)
            Hz.set(i,j,hz)

    Hx.real().show_matrix("%sx_re" % (filename));  Hx.imag().show_matrix("%sx_im" % (filename))    
    Hy.real().show_matrix("%sy_re" % (filename));  Hy.imag().show_matrix("%sy_im" % (filename))    
    Hz.real().show_matrix("%sz_re" % (filename));  Hz.imag().show_matrix("%sz_im" % (filename))  


