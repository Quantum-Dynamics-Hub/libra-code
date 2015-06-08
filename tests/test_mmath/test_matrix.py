import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")


print "\nTest 1: Importing the library and its content"
print "from cygmmath import *"
from cygmmath import *


print "\n Test2: Setting up matrix also test show_matrix() function"
a = MATRIX(3,3)
for i in range(0,3):
    for j in range(0,3):
        a.set(i,j, 0.1+0.1*(i+j) )

for i in range(0,3):
    for j in range(0,3):
        print i,j, a.get(i,j)

a.show_matrix()


print "\n Test3: Properties"
print "trace = ", a.tr()
print "Determinant = ", a.Determinant()  # 0.1*(0.15 - 0.16) - 0.2*(0.1 - 0.12) + 0.3*(0.08 - 0.09) = -0.001 + 0.004 - 0.003 = 0.0


print "\n Test4: JACOBY"
eva = MATRIX(3,3)
eve = MATRIX(3,3)

a.JACOBY_EIGEN(eva,eve, 1e-10)

print "a = \n"
a.show_matrix()

print "eva = \n"
eva.show_matrix()

print "eve = \n"
eve.show_matrix()

print "Checking \n"
print "LHS = M * EVECT\n"
(a*eve).show_matrix()

print "\nRHS = M * EVECT\n"
(eve*eva).show_matrix()


print "\n Test5: Big JACOBY"
for N in [10, 25, 50, 100, 150]:

    a = MATRIX(N,N)
    for i in range(0,N):
        for j in range(0,N):
            a.set(i,j, 0.1+0.1*(i+j) )
    
    eva = MATRIX(N,N)
    eve = MATRIX(N,N)
        
    t = Timer()
    t.start()
    
    a.JACOBY_EIGEN(eva,eve, 1e-10)
    
    t.stop()
    print "Time to diagonalize %5i x %5i matrix = %8.5f sec " % ( N, N, t.show() )
    