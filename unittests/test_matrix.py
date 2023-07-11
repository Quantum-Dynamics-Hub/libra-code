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


class TestMATRIX(unittest.TestCase):
    """ Summary of the tests:
    1 - constructor and members
    2 - set(i,j, val) / get(i,j)
    2a - set(-1, j, val) / get(i,j)
    2b - set(i, -1, val) / get(i,j)
    2c - 2a and 2b together
    2d - set(-1,-1, val) / get(i,j)
    3 - set(i, val) / get(i,j) and get(i)
    4 -  tr / sum
    5 - show_matrix / show_matrix(file)
    6 - bin_dump / bin_load
    7 - col / row 
    8 - T()
    9 - Transpose()
    10 - Assignment and copy constructor
    11 - Diag (versions)
    12 - Init and Init_Unit_Matrix
    13 - swap_cols & swap_rows
    14 - increment/decrement operator, add/scale
    15 - matrix addition and subtraction
    16 - matrix multiplication and division (by a number)
    17 - matrix-matrix multiplication and dot product
    18 - Properties of the matrix
    """

    def test_1(self):
        """Test constructors"""

        print "Checking matrix properties"

        Y = MATRIX()
        X = MATRIX(2,2)

        self.assertEqual(Y.num_of_cols, 0)
        self.assertEqual(Y.num_of_rows, 0)
        self.assertEqual(Y.num_of_elems, 0)

        self.assertEqual(X.num_of_cols, 2)
        self.assertEqual(X.num_of_rows, 2)
        self.assertEqual(X.num_of_elems, 4)



    def test_2(self):        
        """Test set and get"""
        print "Checking set/get functions. Case 1"

        X = MATRIX(2,2)
        X.set(0,0, 1.0); X.set(0,1, -0.5);
        X.set(1,0, 0.6); X.set(1,1,  1.0);

        self.assertEqual( X.get(0,0), 1.0 )
        self.assertEqual( X.get(0,1), -0.5 )
        self.assertEqual( X.get(1,0), 0.6 )
        self.assertEqual( X.get(1,1), 1.0 )

        self.assertEqual( X.get(0), 1.0 )
        self.assertEqual( X.get(1), -0.5 )
        self.assertEqual( X.get(2), 0.6 )
        self.assertEqual( X.get(3), 1.0 )
 

    def test_2a(self):        
        """Test set"""
        print "Checking versions of set: a"

        X = MATRIX(3,3)
        print X.num_of_rows
        print X.num_of_cols
        print X.num_of_elems

        X.set(-1, 1, 1.0); 

        self.assertEqual( X.get(0,1), 1.0 )
        self.assertEqual( X.get(1,1), 1.0 )
        self.assertEqual( X.get(2,1), 1.0 )

        self.assertEqual( X.get(0,0), 0.0 )
        self.assertEqual( X.get(1,0), 0.0 )
        self.assertEqual( X.get(2,0), 0.0 )

        self.assertEqual( X.get(0,2), 0.0 )
        self.assertEqual( X.get(1,2), 0.0 )
        self.assertEqual( X.get(2,2), 0.0 )


    def test_2b(self):        
        """Test set"""
        print "Checking versions of set: b"

        X = MATRIX(3,3)
        X.set(0,-1, 1.4); 

        self.assertEqual( X.get(0,0), 1.4 )
        self.assertEqual( X.get(0,1), 1.4 )
        self.assertEqual( X.get(0,2), 1.4 )

        self.assertEqual( X.get(1,0), 0.0 )
        self.assertEqual( X.get(1,1), 0.0 )
        self.assertEqual( X.get(1,2), 0.0 )

        self.assertEqual( X.get(2,0), 0.0 )
        self.assertEqual( X.get(2,1), 0.0 )
        self.assertEqual( X.get(2,2), 0.0 )


    def test_2c(self):        
        """Test set"""
        print "Checking versions of set: c"

        X = MATRIX(3,3)
        X.set(0,-1, 1.4); 
        X.set(-1,0, -1.0); 

        #  -1.0  1.4  1.4
        #  -1.0  0.0  0.0
        #  -1.0  0.0  0.0


        self.assertEqual( X.get(0,0), -1.0 )
        self.assertEqual( X.get(0,1), 1.4 )
        self.assertEqual( X.get(0,2), 1.4 )

        self.assertEqual( X.get(1,0), -1.0 )
        self.assertEqual( X.get(1,1), 0.0 )
        self.assertEqual( X.get(1,2), 0.0 )

        self.assertEqual( X.get(2,0), -1.0 )
        self.assertEqual( X.get(2,1), 0.0 )
        self.assertEqual( X.get(2,2), 0.0 )


    def test_2d(self):        
        """Test set"""
        print "Checking versions of set"

        X = MATRIX(2,2)
        X.set(-1,-1, 1.5); 

        self.assertEqual( X.get(0,0), 1.5 )
        self.assertEqual( X.get(0,1), 1.5 )

        self.assertEqual( X.get(1,0), 1.5 )
        self.assertEqual( X.get(1,1), 1.5 )



    def test_3(self):        
        """Test set and get"""
        print "Checking set/get functions. Case 2"

        X = MATRIX(2,2)
        X.set(0, 1.0); X.set(1, -0.5);
        X.set(2, 0.6); X.set(3,  1.0);

        self.assertEqual( X.get(0,0), 1.0 )
        self.assertEqual( X.get(0,1), -0.5 )
        self.assertEqual( X.get(1,0), 0.6 )
        self.assertEqual( X.get(1,1), 1.0 )

        self.assertEqual( X.get(0), 1.0 )
        self.assertEqual( X.get(1), -0.5 )
        self.assertEqual( X.get(2), 0.6 )
        self.assertEqual( X.get(3), 1.0 )



    def test_4(self):
        """Test matrix properties"""
        print "Matrix properties"

        X = MATRIX(2,2)
        X.set(0, 1.0); X.set(1, -0.5);
        X.set(2, 0.6); X.set(3,  1.0);

        self.assertEqual( X.tr(), 2.0 )
        self.assertEqual( X.sum(), 2.1 )


    def test_5(self):
        """Test printout"""
        print "Matrix printout"

        X = MATRIX(2,2)
        X.set(0, 1.0); X.set(1, -0.5);
        X.set(2, 0.6); X.set(3,  1.0);

        X.show_matrix()
        X.show_matrix("X.txt")
              
        self.assertEqual( os.path.isfile("X.txt"), True )

        if os.path.isfile("X.txt"):
            f = open("X.txt", "r")
            A = f.readlines()
            f.close()
 
            x = []
            for a in A: 
                tmp = a.split()
                xx = []
                for it in tmp:
                    xx.append( float(it) )
                x.append(xx)

            self.assertEqual( x[0], [1.0, -0.5] )
            self.assertEqual( x[1], [0.6,  1.0] )

            os.system("rm X.txt")


    def test_6(self):
        """Test binary save/load"""
        print "Matrix binary"

        X = MATRIX(2,2)
        X.set(0, 1.0); X.set(1, -0.5);
        X.set(2, 0.6); X.set(3,  1.0);

        X.bin_dump("X.dat")
              
        self.assertEqual( os.path.isfile("X.dat"), True )

        if os.path.isfile("X.dat"):

            Y = MATRIX(2,2)
            Y.bin_load("X.dat")

            self.assertEqual( Y.get(0,0), 1.0 )
            self.assertEqual( Y.get(0,1), -0.5 )
            self.assertEqual( Y.get(1,0), 0.6 )
            self.assertEqual( Y.get(1,1), 1.0 )

            os.system("rm X.dat")


    def test_7(self):
        """Test col and row extraction"""
        print "Col and row"

        X = MATRIX(2,2)
        X.set(0, 1.0); X.set(1, -0.5);
        X.set(2, 0.6); X.set(3,  1.1);

        x1 = X.col(0)
        x2 = X.col(1)  
        y1 = X.row(0)
        y2 = X.row(1)

        self.assertEqual( x1.num_of_cols, 1 )
        self.assertEqual( x1.num_of_rows, 2 )
        self.assertEqual( x1.num_of_elems, 2 )

        self.assertEqual( x2.num_of_cols, 1 )
        self.assertEqual( x2.num_of_rows, 2 )
        self.assertEqual( x2.num_of_elems, 2 )

        self.assertEqual( y1.num_of_cols, 2 )
        self.assertEqual( y1.num_of_rows, 1 )
        self.assertEqual( y1.num_of_elems, 2 )

        self.assertEqual( y2.num_of_cols, 2 )
        self.assertEqual( y2.num_of_rows, 1 )
        self.assertEqual( y2.num_of_elems, 2 )


        self.assertEqual( x1.get(0), 1.0 )
        self.assertEqual( x1.get(1), 0.6 )
        self.assertEqual( x2.get(0), -0.5 )
        self.assertEqual( x2.get(1), 1.1 )

        self.assertEqual( y1.get(0), 1.0 )
        self.assertEqual( y1.get(1), -0.5 )
        self.assertEqual( y2.get(0), 0.6 )
        self.assertEqual( y2.get(1), 1.1 )


        self.assertEqual( type(x1), type(X) )
        self.assertEqual( type(x2), type(X) )


    def test_8(self):
        """Test T()"""
        print "T()"

        X = MATRIX(2,3)
        X.set(0, 1.0); X.set(1, -0.5); X.set(2, 0.5);
        X.set(3, 0.6); X.set(4,  1.1); X.set(5, 0.8);

        Y = X.T()

        self.assertEqual( type(Y), type(X) )

        self.assertEqual( Y.num_of_cols, 2 )
        self.assertEqual( Y.num_of_rows, 3 )
        self.assertEqual( Y.num_of_elems, 6 )

        self.assertEqual( Y.get(0,0), 1.0 )
        self.assertEqual( Y.get(1,0), -0.5 )
        self.assertEqual( Y.get(2,0), 0.5 )

        self.assertEqual( Y.get(0,1), 0.6 )
        self.assertEqual( Y.get(1,1), 1.1 )
        self.assertEqual( Y.get(2,1), 0.8 )


    def test_9(self):
        """Test Transpose()"""
        print "Transpose()"

        X = MATRIX(2,3)
        X.set(0, 1.0); X.set(1, -0.5); X.set(2, 0.5);
        X.set(3, 0.6); X.set(4,  1.1); X.set(5, 0.8);

        X.Transpose()

        self.assertEqual( X.num_of_cols, 2 )
        self.assertEqual( X.num_of_rows, 3 )
        self.assertEqual( X.num_of_elems, 6 )

        self.assertEqual( X.get(0,0), 1.0 )
        self.assertEqual( X.get(1,0), -0.5 )
        self.assertEqual( X.get(2,0), 0.5 )

        self.assertEqual( X.get(0,1), 0.6 )
        self.assertEqual( X.get(1,1), 1.1 )
        self.assertEqual( X.get(2,1), 0.8 )


    def test_10(self):        
        """Test assignment and copy constructor"""
        print "Assignment and copy constructor"

        X = MATRIX(1,2)
        X.set(0, 1.0); X.set(1, -0.5);
        
        Y = X  # Assignment is done by reference!        
        Z = MATRIX(X) # If you want a separate copy - use the copy constructor!
        Z.set(0, 3.0) # this should not change X
        Y.set(0, 2.0) # this will change Y and X, but not Z

        # Check that
        self.assertEqual( Y.get(0,0), 2.0 )
        self.assertEqual( Y.get(0,1), -0.5 )
        self.assertEqual( X.get(0,0), 2.0 )
        self.assertEqual( X.get(0,1), -0.5 )
        self.assertEqual( Z.get(0,0), 3.0 )


    def test_11(self):        
        """Test diag"""
        print "Diag"

        X = MATRIX(3,3)
        X.diag(2, 1.0)

        self.assertEqual( X.get(0,0), 1.0 )
        self.assertEqual( X.get(1,1), 1.0 )
        self.assertEqual( X.get(0,1), 0.0 )
        self.assertEqual( X.get(2,2), 0.0 )
        self.assertEqual( X.get(0,2), 0.0 )
        self.assertEqual( X.get(1,2), 0.0 )

        X = MATRIX(2,3)
        X.diag(3, 2.0)

        self.assertEqual( X.get(0,0), 2.0 )
        self.assertEqual( X.get(1,1), 2.0 )
        self.assertEqual( X.get(1,0), 0.0 )


        X = MATRIX(3,3)
        X.diag(-1.0)

        self.assertEqual( X.get(0,0), -1.0 )
        self.assertEqual( X.get(1,1), -1.0 )
        self.assertEqual( X.get(2,2), -1.0 )


    def test_12(self):        
        """Test Init"""
        print "Init"

        X = MATRIX(2,2)
        X.Init(-1.0)

        self.assertEqual( X.get(0,0), -1.0 )
        self.assertEqual( X.get(0,1), -1.0 )
        self.assertEqual( X.get(1,0), -1.0 )
        self.assertEqual( X.get(1,1), -1.0 )



        X.Init_Unit_Matrix(2.0)

        self.assertEqual( X.get(0,0), 2.0 )
        self.assertEqual( X.get(1,1), 2.0 )
        self.assertEqual( X.get(0,1), 0.0 )
        self.assertEqual( X.get(1,0), 0.0 )


    def test_13(self):
        """swap_cols, swap_rows"""
        print "swap_cols, swap_rows"

        X = MATRIX(3,3)
        X.set(0,0, 1.0);  X.set(0,1, -0.5);  X.set(0,2, -1.5);
        X.set(1,0, 0.6);  X.set(1,1,  1.1);  X.set(1,2,  2.5);
        X.set(2,0, -0.7); X.set(2,1, -1.2);  X.set(2,2,  0.5);

        X1 = MATRIX(X)
        X1.swap_cols(0,1)

        self.assertEqual( X1.get(0,0), -0.5 )
        self.assertEqual( X1.get(1,0), 1.1 )
        self.assertEqual( X1.get(2,0), -1.2 )

        self.assertEqual( X1.get(0,1), 1.0 )
        self.assertEqual( X1.get(1,1), 0.6 )
        self.assertEqual( X1.get(2,1), -0.7 )

        self.assertEqual( X1.get(0,2), -1.5 )
        self.assertEqual( X1.get(1,2), 2.5 )
        self.assertEqual( X1.get(2,2), 0.5 )


        X1.swap_cols(0,2)

        self.assertEqual( X1.get(0,0), -1.5 )
        self.assertEqual( X1.get(1,0), 2.5 )
        self.assertEqual( X1.get(2,0), 0.5 )

        self.assertEqual( X1.get(0,1), 1.0 )
        self.assertEqual( X1.get(1,1), 0.6 )
        self.assertEqual( X1.get(2,1), -0.7 )

        self.assertEqual( X1.get(0,2), -0.5 )
        self.assertEqual( X1.get(1,2), 1.1 )
        self.assertEqual( X1.get(2,2), -1.2 )


        X2 = MATRIX(X)
        X2.swap_rows(0,1)

        self.assertEqual( X2.get(0,0), 0.6 )
        self.assertEqual( X2.get(1,0), 1.0 )
        self.assertEqual( X2.get(2,0), -0.7 )

        self.assertEqual( X2.get(0,1), 1.1 )
        self.assertEqual( X2.get(1,1), -0.5 )
        self.assertEqual( X2.get(2,1), -1.2 )

        self.assertEqual( X2.get(0,2), 2.5 )
        self.assertEqual( X2.get(1,2), -1.5 )
        self.assertEqual( X2.get(2,2), 0.5 )



    def test_14(self):
        """increment/decrement operator, add/scale"""
        print "increment/decrement operator, add/scale"

        X = MATRIX(2,2)
        X.set(0,0, 1.0); X.set(0,1, -0.5); 
        X.set(1,0, 0.6); X.set(1,1,  1.1); 

        X += 5.0

        self.assertAlmostEqual( X.get(0,0), 6.0 )
        self.assertAlmostEqual( X.get(0,1), 4.5 )
        self.assertAlmostEqual( X.get(1,0), 5.6 )
        self.assertAlmostEqual( X.get(1,1), 6.1 )

        X -= 4.5
        self.assertAlmostEqual( X.get(0,0), 1.5 )
        self.assertAlmostEqual( X.get(0,1), 0.0 )
        self.assertAlmostEqual( X.get(1,0), 1.1 )
        self.assertAlmostEqual( X.get(1,1), 1.6 )

        X *= 2.0
        self.assertAlmostEqual( X.get(0,0), 3.0 )
        self.assertAlmostEqual( X.get(0,1), 0.0 )
        self.assertAlmostEqual( X.get(1,0), 2.2 )
        self.assertAlmostEqual( X.get(1,1), 3.2 )

        X /= 3.0
        self.assertAlmostEqual( X.get(0,0), 1.0 )
        self.assertAlmostEqual( X.get(0,1), 0.0 )
        self.assertAlmostEqual( X.get(1,0), 2.2/3.0 )
        self.assertAlmostEqual( X.get(1,1), 3.2/3.0 )

        X.add(1,0, -2.2/3.0)
        self.assertAlmostEqual( X.get(0,0), 1.0 )
        self.assertAlmostEqual( X.get(0,1), 0.0 )
        self.assertAlmostEqual( X.get(1,0), 0.0 )
        self.assertAlmostEqual( X.get(1,1), 3.2/3.0 )

        X.scale(1,1, 3.0/3.2)
        self.assertAlmostEqual( X.get(0,0), 1.0 )
        self.assertAlmostEqual( X.get(0,1), 0.0 )
        self.assertAlmostEqual( X.get(1,0), 0.0 )
        self.assertAlmostEqual( X.get(1,1), 1.0 )


    def test_15(self):
        """matrix addition and subtraction"""
        print "matrix addition and subtraction"

        X = MATRIX(2,2)
        X.set(0,0, 1.0); X.set(0,1, -0.5); 
        X.set(1,0, 0.6); X.set(1,1,  1.1); 

        Y = MATRIX(2,2)
        Y.set(0,0, 0.5); Y.set(0,1, -1.5); 
        Y.set(1,0, 0.2); Y.set(1,1,  1.0); 


        Z = X + Y
        self.assertAlmostEqual( Z.get(0,0), 1.5 )
        self.assertAlmostEqual( Z.get(0,1), -2.0 )
        self.assertAlmostEqual( Z.get(1,0), 0.8 )
        self.assertAlmostEqual( Z.get(1,1), 2.1 )


        Z = X - Y
        self.assertAlmostEqual( Z.get(0,0), 0.5 )
        self.assertAlmostEqual( Z.get(0,1), 1.0 )
        self.assertAlmostEqual( Z.get(1,0), 0.4 )
        self.assertAlmostEqual( Z.get(1,1), 0.1 )

        Z = X + 1.0
        self.assertAlmostEqual( Z.get(0,0), 2.0 )
        self.assertAlmostEqual( Z.get(0,1), 0.5 )
        self.assertAlmostEqual( Z.get(1,0), 1.6 )
        self.assertAlmostEqual( Z.get(1,1), 2.1 )

        Z = Y - 1.0
        self.assertAlmostEqual( Z.get(0,0), -0.5 )
        self.assertAlmostEqual( Z.get(0,1), -2.5 )
        self.assertAlmostEqual( Z.get(1,0), -0.8 )
        self.assertAlmostEqual( Z.get(1,1), 0.0 )





    def test_16(self):
        """matrix multiplication and division"""
        print "matrix multiplication and division"

        X = MATRIX(2,2)
        X.set(0,0, 1.0); X.set(0,1, -0.5); 
        X.set(1,0, 0.6); X.set(1,1,  1.1); 

        Y = MATRIX(2,2)
        Y.set(0,0, 0.5); Y.set(0,1, -1.5); 
        Y.set(1,0, 0.2); Y.set(1,1,  1.0); 


        Z = 2.0 * X 
        self.assertAlmostEqual( Z.get(0,0), 2.0 )
        self.assertAlmostEqual( Z.get(0,1), -1.0 )
        self.assertAlmostEqual( Z.get(1,0), 1.2 )
        self.assertAlmostEqual( Z.get(1,1), 2.2 )

        Z =  X * 2.0
        self.assertAlmostEqual( Z.get(0,0), 2.0 )
        self.assertAlmostEqual( Z.get(0,1), -1.0 )
        self.assertAlmostEqual( Z.get(1,0), 1.2 )
        self.assertAlmostEqual( Z.get(1,1), 2.2 )


        Z = Y / 2.0
        self.assertAlmostEqual( Z.get(0,0), 0.25 )
        self.assertAlmostEqual( Z.get(0,1), -0.75 )
        self.assertAlmostEqual( Z.get(1,0), 0.1 )
        self.assertAlmostEqual( Z.get(1,1), 0.5 )


    def test_17(self):
        """matrix - matrix multiplication"""
        print "matrix - matrix multiplication"

        X = MATRIX(2,2)
        X.set(0,0, 1.0); X.set(0,1, -0.5); 
        X.set(1,0, 0.6); X.set(1,1,  1.1); 

        Y = MATRIX(2,2)
        Y.set(0,0, 0.5); Y.set(0,1, -1.5); 
        Y.set(1,0, 0.2); Y.set(1,1,  1.0); 


        Z = X * Y 
        self.assertAlmostEqual( Z.get(0,0), 0.5 - 0.1 )
        self.assertAlmostEqual( Z.get(0,1), -1.5 - 0.5 )
        self.assertAlmostEqual( Z.get(1,0), 0.3 + 0.22 )
        self.assertAlmostEqual( Z.get(1,1), -0.9 + 1.1 )

        Z.set(-1,-1, 0.0)
        self.assertAlmostEqual( Z.get(0,0), 0.0 )
        self.assertAlmostEqual( Z.get(0,1), 0.0 )
        self.assertAlmostEqual( Z.get(1,0), 0.0 )
        self.assertAlmostEqual( Z.get(1,1), 0.0 )

        Z.product(X, Y)
        self.assertAlmostEqual( Z.get(0,0), 0.5 - 0.1 )
        self.assertAlmostEqual( Z.get(0,1), -1.5 - 0.5 )
        self.assertAlmostEqual( Z.get(1,0), 0.3 + 0.22 )
        self.assertAlmostEqual( Z.get(1,1), -0.9 + 1.1 )

        Z.dot_product(X, Y)
        self.assertAlmostEqual( Z.get(0,0), 0.5 )
        self.assertAlmostEqual( Z.get(0,1), 0.75 )
        self.assertAlmostEqual( Z.get(1,0), 0.12 )
        self.assertAlmostEqual( Z.get(1,1), 1.1 )



    def test_18(self):
        """Properties of the matrix"""
        print "Properties of the matrix"

        X = MATRIX(3,3)
        X.set(0,0, 1.0);  X.set(0,1, -0.5);  X.set(0,2, -1.5);
        X.set(1,0, 0.6);  X.set(1,1,  1.1);  X.set(1,2,  2.5);
        X.set(2,0, -0.7); X.set(2,1, -1.2);  X.set(2,2,  0.5);

        self.assertAlmostEqual( X.max_elt(), 2.5 )

        x = X.max_col_elt(0)
        self.assertAlmostEqual( x[0], 0 )
        self.assertAlmostEqual( x[1], 1.0 )

        x = X.max_col_elt(1)
        self.assertAlmostEqual( x[0], 2 )
        self.assertAlmostEqual( x[1], -1.2 )

        x = X.max_col_elt(2)
        self.assertAlmostEqual( x[0], 1 )
        self.assertAlmostEqual( x[1], 2.5 )


        x = X.max_row_elt(0)
        self.assertAlmostEqual( x[0], 2 )
        self.assertAlmostEqual( x[1], -1.5 )

        x = X.max_row_elt(1)
        self.assertAlmostEqual( x[0], 2 )
        self.assertAlmostEqual( x[1], 2.5 )

        x = X.max_row_elt(2)
        self.assertAlmostEqual( x[0], 1 )
        self.assertAlmostEqual( x[1], -1.2 )

        


if __name__=='__main__':
    unittest.main()




