#*********************************************************************************
#* Copyright (C) 2021 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
 This file tests the correctness of the Pickling support for some Libra datatypes
 Pickling is used in the multiprocessing library
"""

import os
import sys
import math
import time
import unittest
import multiprocessing as mp

from liblibra_core import *



"""

This is a Python implementation of the pickling support - one would only need to 
add the .enable_pickling() line


def matrix_getinitargs(self):
    return (self.num_of_rows, self.num_of_cols )

def matrix_getstate(self):
    res = []
    for i in range(self.num_of_elems):
        res.append( self.get(i) )    
    return tuple(res)

def matrix_setstate(self, res):
    sz = len(res)
    for i in range(sz):
        self.set(i, res[i])

MATRIX.__getinitargs__ = matrix_getinitargs
MATRIX.__getstate__ = matrix_getstate
MATRIX.__setstate__ = matrix_setstate
"""




def summ(ilist, pw):    
    print("Starting summ")
    sz = len(ilist)    
    res = 0
    for i in range(sz):
        res += ilist[i]**pw
    
    print("RES = ", res)
    print("Done with summ")
    return res


def run_the_job(func, _var_pool, nthreads):
    t1 = time.time()
    pool = mp.Pool( nthreads )
    res = pool.starmap( func, _var_pool )
    pool.close()
    pool.join()
    t2 = time.time()
    print(F"Total time {t2 - t1}") 
    return res
    

def test_intList():
    # Pickling of inList type

    x1 = Py2Cpp_int([1,2,3,4])
    x2 = Py2Cpp_int([1,2,3,4])
    var_pool = [ ( x1, 1), ( x2, 2)  ]
    
    res = run_the_job(summ, var_pool, 2)

    print(res)


def summ_matrix(mtx, pw):    
    print("Starting summ_matrix")

    res = 0
    for i in range(mtx.num_of_rows):
        for j in range(mtx.num_of_cols):

            res += mtx.get(i,j)**pw
        
    print("RES = ", res)
    print("Done with summ_matrix")
    return res



def test_MATRIX():
    # Pickling of inList type

    y1 = MATRIX(2,2)
    y1.set(0,0, 1.0);  y1.set(0, 1, 2.0);
    y1.set(1,0, 3.0);  y1.set(1, 1, 4.0);
    y2 = MATRIX(y1)

    var_pool = [ ( y1, 1), ( y2, 2) ]

    # Pickling of MATRIX type
    res = run_the_job(summ_matrix, var_pool, 2)

    print(res)


test_intList()
test_MATRIX()
