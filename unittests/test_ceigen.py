#*********************************************************************************
#* Copyright (C) 2015-2017 Alexey V. Akimov
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
import random

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../_build/src/math_linalg")
sys.path.insert(1,cwd+"/../_build/src/math_meigen")

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    #from cyglibra_core import *
    from cyglinalg import *
    from cygmeigen import *

elif sys.platform=="linux" or sys.platform=="linux2":
    #from liblibra_core import *
    from liblinalg import *
    from libmeigen import *




class TestmEigen(unittest.TestCase):
    """ Summary of the tests:
    1 - determinant
    2 - 
    """

    def test_1(self):
        """Det"""
        print "Det"

        H = MATRIX(2,2)
        H.set(0,0, -0.5);  H.set(0,1, 0.2)
        H.set(1,0, 0.2);   H.set(1,1,  2.0)

        cH = CMATRIX(H)

        self.assertAlmostEqual( det(H), -1.04  )
        self.assertAlmostEqual( det(cH), -1.04+0.0j  )

        self.assertAlmostEqual( FullPivLU_det(H), -1.04  )
        self.assertAlmostEqual( FullPivLU_det(cH), -1.04+0.0j  )


    def test_1a(self):
        """Det"""
        print "Det"

        H = MATRIX(3,3)
        H.set(0,0, -0.5);  H.set(0,1, -0.5);  H.set(0,2, 0.2)
        H.set(1,0, -1.5);  H.set(1,1, -1.5);  H.set(1,2, 2.0)
        H.set(2,0,  1.5);  H.set(2,1,  1.5);  H.set(2,2, 3.0)

        cH = CMATRIX(H)

        self.assertAlmostEqual( det(H), 0.0  )
        self.assertAlmostEqual( det(cH), 0.0+0.0j  )

        self.assertAlmostEqual( FullPivLU_det(H), 0.0  )
        self.assertAlmostEqual( FullPivLU_det(cH), 0.0+0.0j  )


    def test_1b(self):
        """Det"""
        print "Det"
        
        N = 3
        H = MATRIX(N,N)
        for i in xrange(N):
            for j in xrange(N):
                H.set(i,j, 1.0/(1.0+i) + 1.0/(1.0+j) );  

        cH = CMATRIX(H)
        print "Det = ", det(H)
        print "Det = ", FullPivLU_det(H)



    def test_2(self):
        """Solve eigen"""
        H = MATRIX(2,2)
        H.set(0,0, -0.001);  H.set(0,1, 0.001)
        H.set(1,0, 0.001);   H.set(1,1,  0.001)

        S = MATRIX(2,2)   
        S.set(0,0, 1.0);   S.set(0,1, 0.5)
        S.set(1,0, 0.5);   S.set(1,1, 1.0)

        cH = CMATRIX(H)
        cS = CMATRIX(S)
        E = CMATRIX(2,2)
        C = CMATRIX(2,2)

        solve_eigen(H, S, E, C, 0)

        # H*C = S*C*E
        LHS = cH*C
        RHS = cS*C*E
        for i in xrange(4):
            self.assertAlmostEqual( LHS.get(i), RHS.get(i)  )

        # Orthogonality:  C.H() * S * C = I
        csc = C.H() * cS * C
        I = CMATRIX(2,2);  I.identity();
        for i in xrange(4):
            self.assertAlmostEqual( csc.get(i), I.get(i)  )

        #print "H*C = \n"; (cH*C).show_matrix()
        #print "S*C*E = \n"; (cS*C*E).show_matrix()
        #print "C.H() * S * C\n"; (C.H() * cS * C).show_matrix()


    def test_2a(self):
        """Solve eigen"""
        N = 100
        H = MATRIX(N,N)
        S = MATRIX(N,N)

        for i in xrange(N):
            for j in xrange(N):
                H.set(i,j, 1.0/(1.0+i) + 1.0/(1.0+j) );  
                S.set(i,j, math.exp(-0.1*(i-j)**2) );  
         
        cH = CMATRIX(H)
        cS = CMATRIX(S)
        E = CMATRIX(N,N)
        C = CMATRIX(N,N)

        solve_eigen(H, S, E, C, 0)

        # H*C = S*C*E
        LHS = cH*C
        RHS = cS*C*E
        for i in xrange(N*N):
            self.assertAlmostEqual( LHS.get(i), RHS.get(i), 5  )

        # Orthogonality:  C.H() * S * C = I
        csc = C.H() * cS * C
        I = CMATRIX(N,N);  I.identity();
        for i in xrange(N*N):
            self.assertAlmostEqual( csc.get(i), I.get(i), 5  )

        #print "H*C = \n"; (cH*C).show_matrix()
        #print "S*C*E = \n"; (cS*C*E).show_matrix()
        #print "C.H() * S * C\n"; (C.H() * cS * C).show_matrix()


    def test_3(self):
        """Inverse"""
        N = 5
        S = CMATRIX(N,N)

        for i in xrange(N):
            for j in xrange(N):
                S.set(i,j, math.exp(-0.1*(i-j)**2), 0.0 );  
         
        Sinv = CMATRIX(N,N)

        info = FullPivLU_rank_invertible(S)
        print "rank = ", info[0]
        print "is_invertible = ", info[1]

        self.assertEqual( info, [N, 1]  )

        inv_matrix(S, Sinv, -1.0)

        #S.show_matrix()
        #Sinv.show_matrix()
        I = CMATRIX(N,N);  I.identity();

        LHS = Sinv * S
        for i in xrange(N*N):
            self.assertAlmostEqual( LHS.get(i), I.get(i), 7  )


    def test_3a(self):
        """Inverse"""
        N = 25
        S = CMATRIX(N,N)

        for i in xrange(N):
            for j in xrange(N):
                S.set(i,j, math.exp(-0.1*(i-j)**2), 0.0 );  
         
        Sinv = CMATRIX(N,N)

        info = FullPivLU_rank_invertible(S)
        print "rank = ", info[0]
        print "is_invertible = ", info[1]

        self.assertEqual( info, [N, 1]  )

        inv_matrix(S, Sinv, -1.0)

        #S.show_matrix()
        #Sinv.show_matrix()
        I = CMATRIX(N,N);  I.identity();

        LHS = Sinv * S
        for i in xrange(N*N):
            self.assertAlmostEqual( LHS.get(i), I.get(i), 7  )




    def test_5(self):
        """sqrt() and 1/sqrt()"""
        N = 5
        S = CMATRIX(N,N)

        for i in xrange(N):
            for j in xrange(N):
                S.set(i,j, math.exp(-0.1*(i-j)**2), 0.0 );  
         
        S_inv = CMATRIX(N,N)
        S_half_inv = CMATRIX(N,N)
        S_half = CMATRIX(N,N)

        sqrt_matrix(S, S_half, S_half_inv)
        inv_matrix(S, Sinv, -1.0)


        I = CMATRIX(N,N);  I.identity();

        lhs1 = S_half_inv * S_half_inv
        lhs2 = S_half * S_half
        lhs3 = S_half * S_half_inv
        lhs4 = S_half_inv * S_half

        for i in xrange(N*N):
            self.assertAlmostEqual( lhs1.get(i), Sinv.get(i), 7  )
            self.assertAlmostEqual( lhs2.get(i), S.get(i), 7  )
            self.assertAlmostEqual( lhs3.get(i), I.get(i), 7  )
            self.assertAlmostEqual( lhs4.get(i), I.get(i), 7  )






    def test_3b(self):
        """Inverse"""
        N = 100
        S = CMATRIX(N,N)

        for i in xrange(N):
            for j in xrange(N):
                S.set(i,j, math.exp(-0.1*(i-j)**2), 0.0 );  
         
        Sinv = CMATRIX(N,N)

        info = FullPivLU_rank_invertible(S)
        print "rank = ", info[0]
        print "is_invertible = ", info[1]

        self.assertEqual( info, [N, 1]  )

        FullPivLU_inverse(S, Sinv)

        #S.show_matrix()
        #Sinv.show_matrix()
        I = CMATRIX(N,N);  I.identity();

        LHS = Sinv * S
        for i in xrange(N*N):
            self.assertAlmostEqual( LHS.get(i), I.get(i), 5  )



    def test_5(self):
        """Linsys"""
        N = 3 
        M = 2
        C = MATRIX(N,N)
        D = MATRIX(N,M)
        X = MATRIX(N,M)

        for i in xrange(N):
            for j in xrange(N):
                C.set(i,j, random.uniform(-10.0, 10.0) );  
        for i in xrange(M):
            for j in xrange(N):
                D.set(i,j, random.uniform(-1.0, 1.0) );  

        # Solve for X:  C*X = D
        eps = 1e-9
        maxiter = 10000
        linsys_solver(C, D, X, eps, maxiter)
         
        LHS = C * X
 
        for i in xrange(N*M):
            self.assertAlmostEqual( LHS.get(i), D.get(i), 5  )
        


    def no_test_5a(self):
        """Linsys: But this is not very efficient, so don't use this solver for something bigger"""
        N = 30 
        M = 20
        C = MATRIX(N,N)
        D = MATRIX(N,M)
        X = MATRIX(N,M)

        for i in xrange(N):
            for j in xrange(N):
                C.set(i,j, random.uniform(-10.0, 10.0) );  
        for i in xrange(M):
            for j in xrange(N):
                D.set(i,j, random.uniform(-1.0, 1.0) );  

        # Solve for X:  C*X = D
        eps = 1e-9
        maxiter = 10000
        linsys_solver(C, D, X, eps, maxiter)
         
        LHS = C * X
 
        for i in xrange(N*M):
            self.assertAlmostEqual( LHS.get(i), D.get(i), 5  )
        

        


if __name__=='__main__':
    unittest.main()

