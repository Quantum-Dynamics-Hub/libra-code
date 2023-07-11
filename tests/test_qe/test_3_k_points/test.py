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
## \file test.py
# This example demonstrates how to read different k-points and other needed info

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



info, all_e = QE_methods.read_qe_index("x.export/index.xml", [1,2], 1)

c1 = QE_methods.read_qe_wfc("x.export/wfc.1", [1,2], 1)
c2 = QE_methods.read_qe_wfc("x.export/wfc.2", [1,2], 1)
c3 = QE_methods.read_qe_wfc("x.export/wfc.3", [1,2], 1)


s1  = c1.H() * c1
s2  = c2.H() * c2
s3  = c3.H() * c3
print "s1"; s1.show_matrix()
print "s2"; s2.show_matrix()
print "s3"; s3.show_matrix()

g1 = QE_methods.read_qe_wfc_grid("x.export/grid.1", 1)
g2 = QE_methods.read_qe_wfc_grid("x.export/grid.2", 1)
g3 = QE_methods.read_qe_wfc_grid("x.export/grid.3", 1)



### The following operations will not work, because there may be different number of 
### planewaves needed for each k-point. The proper way to compute the overlaps is 
### via additional matrix with integrals over planewaves

#s12  = c1.H() * c2;  print "s12"; s12.show_matrix()
#s13  = c1.H() * c3;  print "s13"; s13.show_matrix()
#s23  = c2.H() * c3;  print "s23"; s23.show_matrix()



print "\n"


