#*********************************************************************************
#* Copyright (C) 2024 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
 This is an extended test suite for validating the reordering algorithm implemented in C++
"""

import os
import sys
import math
import copy
import pytest
#import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

X = MATRIX(1,3); X.set(0,0, 0.0); X.set(0,1, 1.0); X.set(0,2, 2.0);

perm1 = [0,1,2]; X1 = MATRIX(1,3);  X1.set(0,0, 0.0); X1.set(0,1, 1.0); X1.set(0,2, 2.0);
perm2 = [1,0,2]; X2 = MATRIX(1,3);  X2.set(0,0, 1.0); X2.set(0,1, 0.0); X2.set(0,2, 2.0);
perm3 = [2,1,0]; X3 = MATRIX(1,3);  X3.set(0,0, 2.0); X3.set(0,1, 1.0); X3.set(0,2, 0.0);
perm4 = [2,0,1]; X4 = MATRIX(1,3);  X4.set(0,0, 2.0); X4.set(0,1, 0.0); X4.set(0,2, 1.0);
perm5 = [1,2,0]; X5 = MATRIX(1,3);  X5.set(0,0, 1.0); X5.set(0,1, 2.0); X5.set(0,2, 0.0);
#perm6 = [0,1,2]; X6 = MATRIX(1,3);  X6.set(0,0, 2.0); X6.set(0,1, 1.0); X6.set(0,2, 0.0); this is intentionally incorrect

set1 = [ [ perm1, X1], [perm2, X2], [perm3, X3], [perm4, X4], [perm5, X5] ]

@pytest.mark.parametrize(('perm', 'res_mtx'), set1)
def test_permute(perm, res_mtx):

    tX = MATRIX(X)
    tX.permute_cols( Py2Cpp_int(perm) )
    tX.show_matrix();
    res_mtx.show_matrix()
    print(perm)

    assert (tX - res_mtx).get(0,0) < 1e-10 and (tX - res_mtx).get(0,1) < 1e-10 and (tX - res_mtx).get(0,2) < 1e-10


