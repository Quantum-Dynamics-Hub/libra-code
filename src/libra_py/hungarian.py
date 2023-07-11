#*********************************************************************************
#* Copyright (C) 2018  Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: hungarian
   :platform: Unix, Windows
   :synopsis: 
       This module implements the Munkres-Kuhn algorithm for the assignment problem
       The implementation of the code is based on the instructions from:
       http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html

       References:

           1. http://www.public.iastate.edu/~ddoty/HungarianAlgorithm.html
           2. Harold W. Kuhn. The Hungarian Method for the assignment problem.
               *Naval Research Logistics Quarterly*, 2:83-97, 1955.
           3. Harold W. Kuhn. Variants of the Hungarian method for assignment
               problems. *Naval Research Logistics Quarterly*, 3: 253-258, 1956.
           4. Munkres, J. Algorithms for the Assignment and Transportation Problems.
               *Journal of the Society of Industrial and Applied Mathematics*,
               5(1):32-38, March, 1957.
           5. http://en.wikipedia.org/wiki/Hungarian_algorithm

.. moduleauthor:: Alexey V. Akimov

"""


import cmath
import math
import os
import sys

import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
#from libra_py import *



def step1(X):
    """
    For each row of the matrix, find the smallest element and subtract it from every 
    element in its row. Go to Step 2. 
    """

    n = X.num_of_cols

    for i in range(0,n):
        res = X.min_row_elt(i)
        for j in range(0,n):
            X.add(i, j, -res[1])  

    return 2


def step2(X, M, Rcov, Ccov):
    """
    Find a zero (Z) in the resulting matrix. If there is no starred zero in its row or column, 
    star Z. Repeat for each element in the matrix. Go to Step 3. 
    """

    n = X.num_of_cols

    for i in range(0,n):
        for j in range(0,n):
            if math.fabs(X.get(i,j))<1e-10 and Rcov[i]==0 and Ccov[j]==0:
                M[i][j] = 1
                Rcov[i] = 1
                Ccov[j] = 1

    """ 
    Before we go on to Step 3, we uncover all rows and columns so that we can use the
    cover vectors to help us count the number of starred zeros
    """
    for i in range(0,n):
        Rcov[i] = 0
        Ccov[i] = 0
    
    return 3
    

def step3(M, Ccov):
    """   
    Cover each column containing a starred zero.  
    If K columns are covered, the starred zeros describe a complete set of unique assignments.
    In this case, Go to DONE, otherwise, Go to Step 4.
    """

    n = len(M)

    for i in range(0,n):
        for j in range(0,n):
            if M[i][j]==1:
                Ccov[j] = 1

    colcount = 0
    for i in range(0,n):
        colcount += Ccov[i]

    step = 4
    if colcount>=n:
        step = 7

    return step, colcount



def find_a_zero(X, Rcov, Ccov):

    n = X.num_of_cols
    done = False
    row = -1
    col = -1

    r = 0
    while(not done):

        c = 0
        while(True):
            if(math.fabs(X.get(r,c))<1e-20 and Rcov[r]==0 and Ccov[c]==0):
                row, col, done = r, c, True

            c += 1
            if(c>=n or done):
                break
         
        r += 1
        if(r>=n):
            done = True

    return row, col    


def find_star_in_row(M, row):
    n = len(M)

    col = -1
    for i in range(0,n):
        if M[row][i] == 1:
            col = i

    return col

def find_star_in_col(M, col):
    n = len(M)

    row = -1
    for i in range(0,n):
        if M[i][col] == 1:
            row = i

    return row


def find_prime_in_row(M, row):
    n = len(M)

    col = -1
    for i in range(0,n):
        if M[row][i] == 2:
            col = i

    return col




def step4(X, M, Rcov, Ccov):
    """
    Find a noncovered zero and prime it. 
    If there is no starred zero in the row containing this primed zero, Go to Step 5. 
    Otherwise, cover this row and uncover the column containing the starred zero. 
    Continue in this manner until there are no uncovered zeros left. 
    Save the smallest uncovered value and Go to Step 6.
    """

    done = False
    row = -1
    col = -1
    step = 6
    path_row_0, path_col_0 = None, None

    while(not done):    
        row, col = find_a_zero(X,Rcov, Ccov)
        if row==-1:
            done = True
            step = 6
        else:
            M[row][col] = 2

            col_ = find_star_in_row(M, row)

            if col_>-1:
                col = col_
                Rcov[row] = 1
                Ccov[col] = 0
            else:
                done = True
                path_row_0 = row
                path_col_0 = col
                step = 5

    return step, path_row_0, path_col_0


def augment_path(M, path):

    for p in path:
        if M[p[0]][p[1]]==1:
            M[p[0]][p[1]] = 0
        else:
            M[p[0]][p[1]]=1


def clear_covers(Rcov, Ccov):

    n = len(Rcov)
    for i in range(0,n):
        Rcov[i] = 0
        Ccov[i] = 0


def erase_primes(M):

    n = len(M)

    for i in range(0,n):
        for j in range(0,n):
            if M[i][j]==2:
                M[i][j] = 0


def step5(M, path_row_0, path_col_0, Rcov, Ccov):
    """
    Construct a series of alternating primed and starred zeros as follows.
    Let Z0 represent the uncovered primed zero found in Step 4.
    Let Z1 denote the starred zero in the column of Z0 (if any). 
    Let Z2 denote the primed zero in the row of Z1 (there will always be one).  
    Continue until the series terminates at a primed zero that has no starred zero in its column.  
    Unstar each starred zero of the series, star each primed zero of the series, erase all primes
    and uncover every line in the matrix. Return to Step 3. 
    """ 
    done = False
    r = -1
    c = -1

    path = [[path_row_0, path_col_0]]
    path_count = 1

    while(not done):

        #print path
      
        r = find_star_in_col(M, path[path_count-1][1])
        if(r>-1):
            path_count += 1
            path.append([ r, path[path_count-2][1] ])      
        else:
            done = True
      
        if(not done):
            c = find_prime_in_row(M, path[path_count-1][0])
            path_count += 1
            path.append([ path[path_count-2][0], c])      


    augment_path(M, path)
    clear_covers(Rcov, Ccov)  
    erase_primes(M)

    return 3, path



def find_smallest(X, Rcov, Ccov):

    n = X.num_of_cols
    minval = 1.0e+25

    for i in range(0,n):
        for j in range(0,n):
            if(Rcov[i]==0 and Ccov[j]==0):
                if minval>X.get(i,j):
                    minval = X.get(i,j)

    return minval



def step6(X, Rcov, Ccov):
    """
    Add the value found in Step 4 to every element of each covered row, and subtract it from
    every element of each uncovered column.  Return to Step 4 without altering any stars, 
    primes, or covered lines. 
    """

    minval = find_smallest(X, Rcov, Ccov)

    n = X.num_of_cols

    for i in range(0,n):
        for j in range(0,n):
    
            if Rcov[i]==1:
                X.add(i,j, minval)

            if Ccov[j]==0:
                X.add(i,j, -minval)
    return 4


def compute_assignment(M):

    res = []
    n = len(M)

    for i in range(0,n):
        for j in range(0,n):
            if M[i][j] == 1:
                res.append([i,j])
    return res

def show_M(M):

    n = len(M)

    for i in range(0,n):
        res = ""
        for j in range(0,n):
            res = res + " %3i " % M[i][j]
        print(res)
 

def minimize(_X, verbosity=0):

    X = MATRIX(_X)

    n = X.num_of_cols

    M = intList2()
    Rcov = intList()
    Ccov = intList()

    for i in range(0,n):
        Rcov.append(0)
        Ccov.append(0)

        Mi = intList() 
        for j in range(0,n):
            Mi.append(0)
        M.append(Mi)


    done = False
    step = 1
    path_row_0, path_col_0 = None, None
    path = []

    iter_ = 2
    res = []
    while(not done):

        if verbosity > 0:
            print("** Iteration %3i **" % (iter_))
        iter_ += 1

        if step==1:
            if verbosity > 0:
                print("Step 1")
            step = step1(X)

        elif step==2:
            if verbosity > 0:
                print("Step 2")
            step = step2(X, M, Rcov, Ccov)

        elif step==3:
            if verbosity > 0:
                print("Step 3")
            step, colcount = step3(M, Ccov)

        elif step==4:
            if verbosity > 0:
                print("Step 4")
            step, path_row_0, path_col_0 = step4(X, M, Rcov, Ccov)

        elif step==5:
            if verbosity > 0:
                print("Step 5")
            step, path = step5(M, path_row_0, path_col_0, Rcov, Ccov)

        elif step==6:
            if verbosity > 0:
                print("Step 6")
            step = step6(X, Rcov, Ccov)

        elif step==7:
            res = compute_assignment(M)
            if verbosity > 0:
                print("Done")
                print(res)
            done = True

        if verbosity > 0:
            print("X = "); X.show_matrix()    
            print("M = "); show_M(M)
            print("Ccov = ", Cpp2Py(Ccov))
            print("Rcov = ", Cpp2Py(Rcov))
            print("============================")

    return res    


def maximize(_X, verbosity=0):
    """
    Minimize the negative of the original matrix
 
    We also need to shift the negative matrix up rigidly, so all its 
    elements are > 0.0
    """

    X = MATRIX(_X)
    n = X.num_of_cols
    maxval = 0.0
    for i in range(0,n):
        for j in range(0,n):
            if X.get(i,j)>maxval:
                maxval = X.get(i,j)

    for i in range(0,n):
        for j in range(0,n):
            X.scale(i,j, -1.0)
            X.add(i,j, maxval + 1e-5)

    return minimize(X, verbosity)


def _test_setup():

    S = MATRIX(3,3)
    S.set(0,0, 1.0); S.set(0,1, 2.0); S.set(0,2, 3.0); 
    S.set(1,0, 2.0); S.set(1,1, 4.0); S.set(1,2, 6.0); 
    S.set(2,0, 3.0); S.set(2,1, 6.0); S.set(2,2, 9.0); 

    X = MATRIX(4,4)
    X.set(0,0, 1.0); X.set(0,1, 2.0); X.set(0,2, 3.0); X.set(0,3, 4.0); 
    X.set(1,0, 2.0); X.set(1,1, 4.0); X.set(1,2, 6.0); X.set(1,3, 8.0); 
    X.set(2,0, 3.0); X.set(2,1, 6.0); X.set(2,2, 9.0); X.set(2,3,12.0); 
    X.set(3,0, 4.0); X.set(3,1, 8.0); X.set(3,2,12.0); X.set(3,3,16.0); 

    # identical
    # | 1 0 0 0 |
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    a = MATRIX(4,4)
    a.set(0,0,1.0);  a.set(1,1,1.0); 
    a.set(2,2,1.0);  a.set(3,3,1.0); 

    # non-identical (4x4) matrix
    # | 1 0 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    # | 0 1 0 0 |
    b = MATRIX(4,4)
    b.set(0,0,1.0);  b.set(1,2,1.0); 
    b.set(2,3,1.0);  b.set(3,1,1.0);

    # non-identical (4x4) matrix
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 1 0 0 0 |
    # | 0 0 0 1 |
    c = MATRIX(4,4)
    c.set(0,1,1.0);  c.set(1,2,1.0); 
    c.set(2,0,1.0);  c.set(3,3,1.0);

    # non-identical (8x8) matrix
    # | 1 0 0 0 0 0 0 0 |
    # | 0 1 0 0 0 0 0 0 |
    # | 0 0 0 1 0 0 0 0 |
    # | 0 0 1 0 0 0 0 0 |
    # | 0 0 0 0 1 0 0 0 |
    # | 0 0 0 0 0 0 0 1 |
    # | 0 0 0 0 0 1 0 0 |
    # | 0 0 0 0 0 0 1 0 |

    d = MATRIX(8,8)
    d.set(0,0,1.0);  d.set(1,1,1.0);
    d.set(2,3,1.0);  d.set(3,2,1.0);
    d.set(4,4,1.0);  d.set(5,7,1.0);
    d.set(6,5,1.0);  d.set(7,6,1.0);

    return S, X, a, b, c, d



class TestHungarian(unittest.TestCase):

    def test_minimize(self):
        """Tests the reordering algorithm"""
        S, X, a,b,c,d = _test_setup()

        P = minimize(S)
        print("Input matrix "); S.show_matrix()
        print("Assignment map = ", P)
        for p in P:
            self.assertIn(p, [[0,2],[1,1],[2,0]])

        P = minimize(X)
        print("Input matrix "); X.show_matrix()
        print("Assignment map = ", P)
        for p in P:
            self.assertIn(p, [[0,3],[1,2],[2,1],[3,0]])


    def test_maximize(self):
        """Tests the reordering algorithm"""
        S, X, a,b,c,d = _test_setup()

        P = maximize(a)
        print("Input matrix "); a.show_matrix()
        print("Assignment map = ", P)
        for p in P:
            self.assertIn(p, [[0,0],[1,1],[2,2],[3,3]])

        P = maximize(b)
        print("Input matrix "); b.show_matrix()
        print("Assignment map = ", P)
        for p in P:
            self.assertIn(p, [[0,0],[1,2],[2,3],[3,1]])

        P = maximize(c)
        print("Input matrix "); c.show_matrix()
        print("Assignment map = ", P)
        for p in P:
            self.assertIn(p, [[0,1],[1,2],[2,0],[3,3]])

        P = maximize(d)
        print("Input matrix "); d.show_matrix()
        print("Assignment map = ", P)
        for p in P:
            self.assertIn(p, [[0,0],[1,1],[2,3],[3,2],[4,4],[5,7],[6,5],[7,6]])


if __name__=='__main__':
    unittest.main()

