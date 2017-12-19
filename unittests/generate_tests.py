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
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


import PyQuante

# This file will compare the values of the integrals computed with Libra and PyQuante
# and then print out the set of input/output pairs that will be then fixed and used in
# the unit tests. This file is not a unit test - it is just shows how the integrals
# computed by the two packages are related to each other and that they agree, at least
# at the time when this file is created

def run_test(a1,R1, lst1, a2, R2, lst2, tol, mode, f, tests):
# mode - "test" - compare with PyQuante; everythig else - generate the reference data
# f - the file handler
# tests = [0 - unnormalized overlaps, 1 - normalized overlaps, 2 - kinetic energy]

    for n1 in lst1:
        for m1 in lst1:
            for k1 in lst1:
  
                for n2 in lst2:
                    for m2 in lst2:
                        for k2 in lst2:

                            if mode=="test":

                                if 0 in tests:
                                    val_ref = PyQuante.cints.overlap(a1, (n1,m1,k1), (R1.x,R1.y,R1.z),  a2, (n2,m2,k2), (R2.x,R2.y,R2.z) )
  
                                    g1 = PrimitiveG()
                                    g2 = PrimitiveG()
                                    g1.init(n1,m1,k1, a1, R1)
                                    g2.init(n2,m2,k2, a2, R2)

                                    val = gaussian_overlap(g1,g2,0)

                                    if abs(val_ref-val)>tol:
                                        print n1,m1,k1, n2,m2,k2, val_ref, val


                                if 1 in tests:
 
                                    p1 = PyQuante.PGBF.PGBF(a1,(R1.x,R1.y,R1.z),(n1,m1,k1))
                                    p2 = PyQuante.PGBF.PGBF(a2,(R2.x,R2.y,R2.z),(n2,m2,k2))

                                    val_ref = p1.overlap(p2)


                                    g1 = PrimitiveG()
                                    g2 = PrimitiveG()
                                    g1.init(n1,m1,k1, a1, R1)
                                    g2.init(n2,m2,k2, a2, R2)

                                    val = gaussian_overlap(g1,g2)

                                    if abs(val_ref-val)>tol:
                                        print n1,m1,k1, n2,m2,k2, val_ref, val


                                if 2 in tests:

                                    p1 = PyQuante.PGBF.PGBF(a1,(R1.x,R1.y,R1.z),(n1,m1,k1))
                                    p2 = PyQuante.PGBF.PGBF(a2,(R2.x,R2.y,R2.z),(n2,m2,k2))

                                    val_ref = p1.kinetic(p2)


                                    g1 = PrimitiveG()
                                    g2 = PrimitiveG()
                                    g1.init(n1,m1,k1, a1, R1)
                                    g2.init(n2,m2,k2, a2, R2)

                                    val = kinetic_integral(g1,g2)

                                    if abs(val_ref-val)>tol:
                                        print n1,m1,k1, n2,m2,k2, val_ref, val



                            else:

                                g1 = PrimitiveG()
                                g2 = PrimitiveG()
                                g1.init(n1,m1,k1, a1, R1)
                                g2.init(n2,m2,k2, a2, R2)

                                unnorm_ovlp = gaussian_overlap(g1,g2,0)
                                norm_ovlp = gaussian_overlap(g1,g2)
                                kin = kinetic_integral(g1,g2)

                                line =  "%8.5f %8.5f %8.5f %8.5f %i %i %i %8.5f %8.5f %8.5f %8.5f %i %i %i %12.8f %12.8f %12.8f \n" % (a1, R1.x, R1.y, R1.z, n1, m1, k1, a2, R2.x, R2.y, R2.z, n2, m2, k2, unnorm_ovlp, norm_ovlp, kin)
                                f.write(line)



def run_test2(a1, lst1, a2, lst2, tol, mode, f, tests):

    run_test(a1, VECTOR(0.0, 0.0, 0.0), lst1,  a2, VECTOR(1.0, 0.0, 0.0), lst2, tol, mode, f, tests)     
    run_test(a1, VECTOR(0.0, 0.0, 0.0), lst1,  a2, VECTOR(0.0, 1.0, 0.0), lst2, tol, mode, f, tests)
    run_test(a1, VECTOR(0.0, 0.0, 0.0), lst1,  a2, VECTOR(0.0, 0.0, 1.0), lst2, tol, mode, f, tests)
    run_test(a1, VECTOR(0.0, 0.0, 0.0), lst1,  a2, VECTOR(1.0, 0.0, 1.0), lst2, tol, mode, f, tests)
    run_test(a1, VECTOR(0.0, 0.0, 0.0), lst1,  a2, VECTOR(0.0, 1.0, 1.0), lst2, tol, mode, f, tests)
    run_test(a1, VECTOR(0.0, 0.0, 0.0), lst1,  a2, VECTOR(1.0, 1.0, 0.0), lst2, tol, mode, f, tests)
    run_test(a1, VECTOR(0.0, 0.0, 0.0), lst1,  a2, VECTOR(1.0, 1.0, 1.0), lst2, tol, mode, f, tests)


def run_test3(lst1, lst2, tol, mode, f, tests):

    run_test2(0.5, lst1, 0.5, lst2, tol, mode, f, tests) 
    run_test2(0.5, lst1, 1.0, lst2, tol, mode, f, tests)
    run_test2(1.0, lst1, 0.5, lst2, tol, mode, f, tests)


# Unnormalized overlaps
run_test3([0,1,2,3,4], [0,1,2,3,4], 1e-8, "test", None, [0])

# Normalized overlaps
run_test3([0,1,2,3,4], [0,1,2,3,4], 1e-8, "test", None, [1])

# Kinetic integrals
run_test3([0,1,2], [0,1,2], 1e-8, "test", None, [2])


# Reference integrals
f = open("reference_integrals.txt", "w")
run_test3([0,1,2,3], [0,1,2,3], 1e-8, "no-test", f, [0,1,2])
f.close()






