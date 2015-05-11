import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
#sys.path.insert(1,cwd+"/../../_build/src/cell")
sys.path.insert(1,cwd+"/../../_build/src/hamiltonian")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cyghamiltonian import *

print "\nTest 2: SAC_Ham"
f = open("sac_ham.txt","w")
f1 = open("sac_ham1.txt","w")
f2 = open("sac_ham2.txt","w")
for i in range(-100,100):
    x = 0.1 * i
    res = SAC_Ham(x,[])  # same as using [0.010, 1.600, 0.005, 1.00]
    H = res[1]
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H.get(0,0), H.get(1,1), H.get(0,1) ) )
    H1 = res[2]
    f1.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H1.get(0,0), H1.get(1,1), H1.get(0,1) ) )
    H2 = res[3]
    f2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H2.get(0,0), H2.get(1,1), H2.get(0,1) ) )

f.close()
f1.close()
f2.close()


print "\nTest 2.1: SAC_Ham with alternative parameters"
f = open("m-sac_ham.txt","w")
f1 = open("m-sac_ham1.txt","w")
f2 = open("m-sac_ham2.txt","w")
for i in range(-100,100):
    x = 0.1 * i
    res = SAC_Ham(x,[0.10, 1.600, 0.005, 1.00])  # default: [0.010, 1.600, 0.005, 1.00]
    H = res[1]
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H.get(0,0), H.get(1,1), H.get(0,1) ) )
    H1 = res[2]
    f1.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H1.get(0,0), H1.get(1,1), H1.get(0,1) ) )
    H2 = res[3]
    f2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H2.get(0,0), H2.get(1,1), H2.get(0,1) ) )

f.close()
f1.close()
f2.close()




print "\nTest 3: DAC_Ham"
f = open("dac_ham.txt","w")
f1 = open("dac_ham1.txt","w")
f2 = open("dac_ham2.txt","w")
for i in range(-100,100):
    x = 0.1 * i
    res = DAC_Ham(x,[])
    H = res[1]
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H.get(0,0), H.get(1,1), H.get(0,1) ) )
    H1 = res[2]
    f1.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H1.get(0,0), H1.get(1,1), H1.get(0,1) ) )
    H2 = res[3]
    f2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H2.get(0,0), H2.get(1,1), H2.get(0,1) ) )

f.close()
f1.close()
f2.close()



print "\nTest 4: ECWR_Ham"
f = open("ecwr_ham.txt","w")
f1 = open("ecwr_ham1.txt","w")
f2 = open("ecwr_ham2.txt","w")
for i in range(-100,100):
    x = 0.1 * i
    res = ECWR_Ham(x,[])
    H = res[1]
    f.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H.get(0,0), H.get(1,1), H.get(0,1) ) )
    H1 = res[2]
    f1.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H1.get(0,0), H1.get(1,1), H1.get(0,1) ) )
    H2 = res[3]
    f2.write("%8.5f  %8.5f  %8.5f  %8.5f\n" % (res[0], H2.get(0,0), H2.get(1,1), H2.get(0,1) ) )

f.close()
f1.close()
f2.close()



