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


print "\nTest 2: Using timer"
print "t = Timer()"
print "t.start()"
print "x = 0.0"
print "for i in range(0,1000000):"
print "    x = x + math.sin(i*math.pi)"
print "t.stop()"
print "t.show()"

t = Timer()
t.start()
x = 0.0
for i in range(0,1000000):
    x = x + math.sin(i*math.pi)
t.stop()
print "Time to compute = ", t.show(), " sec"


print "t = Timer()"
print "t.start()"
print "x = 0.0"
print "for i in range(0,10000000):"
print "    x = x + math.sin(i*math.pi)"
print "t.stop()"
print "t.show()"
t = Timer()
t.start()
x = 0.0
for i in range(0,10000000):
    x = x + math.sin(i*math.pi)
t.stop()
print "Time to compute = ", t.show(), " sec"

