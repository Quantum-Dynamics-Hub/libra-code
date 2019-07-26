#*********************************************************************************
#* Copyright (C) 2019 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""

 This example demonstrates the functions:

   vector<vector<int> > compute_mapping(vector<vector<int> >& inp, vector<int>& npts);
   int compute_imapping(vector<int>& inp, vector<int>& npts);

   defined in:   dyn/wfcgrid/Grid_functions.h


"""

import os
import sys
import math
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

"""
 Imagine a 3D grid with:
 3 points in the 1-st dimension
 2 points in the 2-nd dimension  
 4 points in the 3-rd dimension
"""

inp = intList2()
npts = Py2Cpp_int([3,2,4])

res = compute_mapping(inp, npts);

print "The number of points = ", len(res)
print "The number of dimensions = ", len(res[0])

cnt = 0
for i in res:
    print "point # ",cnt, Cpp2Py(i)
    print "index of that point in the global array =", compute_imapping(i, Py2Cpp_int([3,2,4]))
    cnt +=1
    
