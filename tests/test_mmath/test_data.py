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

print "\nTest 2: Constructing"
print "x = DATA([1.0, 0.5, 2.0, -0.5])"
x = DATA([1.0, 0.5, 2.0, -0.5])

print "\nTest 3: Estimators"
ave, var, sd, se, mse, mae, rmse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
#print "x.Calculate_Estimators(ave, var, sd, se, mse, mae, rmse)"
#x.Calculate_Estimators(ave, var, sd, se, mse, mae, rmse)
print "ave = ",ave
print "var = ",var
print "sd = ",sd
print "se = ",se
print "mse = ",mse
print "mae = ",mae
print "rmse = ",rmse


