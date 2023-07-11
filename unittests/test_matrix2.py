#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
import cmath
import math
import os
import sys
import unittest

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cyglinalg import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from liblinalg import *



def makeAB():

    N = 5
    A = MATRIX(N,N)
    for i in range(0,N):
        for j in range(0,N):
            A.set(i,j, -2.0*i+j*j )

    return A


class TestMATRIX2(unittest.TestCase):
    """ Summary of the tests:
    1 - extracting and inserting square matrices
    2 - -- general rectangular matrices
    """

    def test_1(self):
        """Extracting/inserting square matrices"""
        A = makeAB()
        
        stenc_list = [ [0,1], [3,4], [1,3], [0,1,2], [0,3,4] ]  

        for stenc in stenc_list:

            print "Extracting with stencil ", stenc
       
            sz = len(stenc)
            a = MATRIX(sz,sz)         
            pop_submatrix(A, a, stenc); 

            B = MATRIX(5,5)
            push_submatrix(B, a, stenc); 

            for i in xrange(sz):
                for j in xrange(sz):

                    self.assertAlmostEqual( a.get(i,j), A.get(stenc[i],stenc[j]) )
                    self.assertAlmostEqual( a.get(i,j), B.get(stenc[i],stenc[j]) )

            print "Extracted matrix"; a.show_matrix() 
            print "Target matrix after insersion"; B.show_matrix()


       
    def test_2(self):
        """Extracting/inserting rectangular matrices"""

        A = makeAB()
        
        stenc_pair_list = [ [[2,4] , [1,3,4]],  [[0,1] , [2,3,4]], [[2,3,4] , [0,1]] ]  

        for stenc in stenc_pair_list:

            print "Extracting with stencil ", stenc[0], stenc[0]
       
            sz1 = len(stenc[0])
            sz2 = len(stenc[1])

            a = MATRIX(sz1,sz2)         
            pop_submatrix(A, a, stenc[0], stenc[1]); 

            B = MATRIX(5,5)
            push_submatrix(B, a, stenc[0], stenc[1]); 

            for i in xrange(sz1):
                for j in xrange(sz2):

                    self.assertAlmostEqual( a.get(i,j), A.get(stenc[0][i],stenc[1][j]) )
                    self.assertAlmostEqual( a.get(i,j), B.get(stenc[0][i],stenc[1][j]) )

            print "Extracted matrix"; a.show_matrix() 
            print "Target matrix after insersion"; B.show_matrix()





if __name__=='__main__':
    unittest.main()




