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


print "\nTest 2: Creating SPACE_GROUP object"
sg = SPACE_GROUP("P_1")
#sg = SPACE_GROUP("Cm_m_2")

sz = len(sg.operators)
print "Number of operators = ", sz

for i in range(0,sz):
    sg.operators[i].show_matrix()


print "\nTest 3: Adding symmetry-equivalent positions: P_1 symmetry"
all_r = VECTORList()
Apply_Symmetry("P_1",VECTOR(0.0, 0.0, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-equivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z


print "\nTest 3a: Adding symmetry-equivalent positions: P_1 symmetry"
all_r = VECTORList()
Apply_Symmetry("P_1",VECTOR(0.1, 0.0, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-equivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z


print "\nTest 3b: Adding symmetry-equivalent positions: P_1 symmetry"
all_r = VECTORList()
Apply_Symmetry("P_1",VECTOR(0.1, 0.21, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-equivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z




print "\nTest 4: Adding symmetry-equivalent positions: Cm_m_2 symmetry"
all_r = VECTORList()
Apply_Symmetry("Cm_m_2",VECTOR(0.0, 0.0, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-equivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z


print "\nTest 4a: Adding symmetry-equivalent positions: Cm_m_2 symmetry"
all_r = VECTORList()
Apply_Symmetry("Cm_m_2",VECTOR(0.1, 0.0, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-equivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z


print "\nTest 4b: Adding symmetry-equivalent positions: Cm_m_2 symmetry"
all_r = VECTORList()
Apply_Symmetry("Cm_m_2",VECTOR(0.1, 0.21, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-equivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z




print "\nTest 5: Adding symmetry-equivalent positions: I_m_-3_m symmetry"
all_r = VECTORList()
Apply_Symmetry("I_m_-3_m",VECTOR(0.0, 0.0, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-equivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z


print "\nTest 5a: Adding symmetry-equivalent positions: I_m_-3_m symmetry"
all_r = VECTORList()
Apply_Symmetry("I_m_-3_m",VECTOR(0.1, 0.0, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-inequivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z


print "\nTest 5b: Adding symmetry-equivalent positions: I_m_-3_m symmetry"
all_r = VECTORList()
Apply_Symmetry("I_m_-3_m",VECTOR(0.1, 0.21, 0.0), all_r)

sz = len(all_r)
print "Number of symmetry-inequivalent positions = ", sz
for i in range(0,sz):
    print i,all_r[i].x, all_r[i].y, all_r[i].z

