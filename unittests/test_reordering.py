#*********************************************************************************  
#* Copyright (C) 2017-2018 Kosuke Sato, Alexey V. Akimov 
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
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
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def get_U(N):
    U = CMATRIX(N,N)
    for i in xrange(N):
        U.set(i,i, 1.0+0.0j) 
    return U

def permute(U, perm):
    pU = CMATRIX(U)
    pU.permute_cols(Py2Cpp_int(perm))
    return pU

def run_test(N, perm):
    U  = get_U(N)
    Ut = permute(U, perm)
    X = U.H() * Ut
    res = Cpp2Py(get_reordering(X))

    print "Ut = "; Ut.show_matrix()
    print "Permutation", res
    print "X = "; X.show_matrix()         

    return res

    

class TestReordering(unittest.TestCase):
    def test_reordering(self):
        """Tests the reordering algorithm"""

        perm = [0,1]    
        res = run_test(2, perm)
        self.assertEqual(res, perm)

        perm = [1,0]    
        res = run_test(2, perm)
        self.assertEqual(res, perm)


        perm = [0,1,2]    
        res = run_test(3, perm)
        self.assertEqual(res, perm)

        perm = [0,2,1]    
        res = run_test(3, perm)
        self.assertEqual(res, perm)

        perm = [1,0,2]    
        res = run_test(3, perm)
        self.assertEqual(res, perm)

        perm = [1,2,0]    
        res = run_test(3, perm)
        self.assertEqual(res, perm)

        perm = [2,1,0]    
        res = run_test(3, perm)
        self.assertEqual(res, perm)

        perm = [2,0,1]    
        res = run_test(3, perm)
        self.assertEqual(res, perm)




        perm = [0,1,2,3]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [0,2,1,3]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [1,0,2,3]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [1,2,0,3]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [2,1,0,3]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [2,0,1,3]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)


        perm = [1,2,3,0]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [1,3,2,0]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [2,1,3,0]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [2,3,1,0]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [3,2,1,0]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [3,1,2,0]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)


        perm = [2,3,0,1]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [2,0,3,1]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [3,2,0,1]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [3,0,2,1]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [0,3,2,1]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [0,2,3,1]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)


        perm = [3,0,1,2]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [3,1,0,2]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [0,3,1,2]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [0,1,3,2]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [1,0,3,2]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)
        perm = [1,3,0,2]    
        res = run_test(4, perm)
        self.assertEqual(res, perm)





if __name__=='__main__':
    unittest.main()

