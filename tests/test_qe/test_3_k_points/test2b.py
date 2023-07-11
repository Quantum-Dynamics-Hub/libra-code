#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file test2b.py
# Here, we will implement the functions for computing overlaps
# These are Python implementation of the functions to compute the overlap integrals with multiple k-points
# This is much faster but still takes some time
# Compose the overall overlap matrix!
# 

import os
import sys
import math
import cmath

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



info, all_e = QE_methods.read_qe_index("x.export/index.xml", [1,2], 1)

print "The total # of k-points is: ", info["nk"]

coeff = []
grid = []

for ik in xrange(info["nk"]):
    print ik, info["k"][ik]
    coeff.append(QE_methods.read_qe_wfc("x.export/wfc.%i" % (ik+1), [1,2], 0))
    grid.append( QE_methods.read_qe_wfc_grid("x.export/grid.%i" % (ik+1) , 0) )


S = CMATRIX(info["nk"]*2, info["nk"]*2)

for ik1 in xrange(info["nk"]):
    for ik2 in range(ik1, info["nk"]):

        s = pw_overlap(info["k"][ik1], info["k"][ik2], coeff[ik1], coeff[ik2], grid[ik1], grid[ik2])
      
        if ik2!=ik1:
            push_submatrix(S, s, range(ik1*2, (ik1+1)*2), range(ik2*2, (ik2+1)*2))
            push_submatrix(S, s.H(), range(ik2*2, (ik2+1)*2), range(ik1*2, (ik1+1)*2))

            print ik1, ik2;  s.show_matrix()
            print ik2, ik1;  s.H().show_matrix()

        else:
            push_submatrix(S, s, range(ik1*2, (ik1+1)*2), range(ik2*2, (ik2+1)*2))

            print ik1, ik2;  s.show_matrix()
        
print "The overall overlap matrix: "; S.show_matrix()



