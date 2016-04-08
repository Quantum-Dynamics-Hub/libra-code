#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
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
cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/mmath")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/util")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygconverters import *
from cygutil import *

print "\nTest2: void show_vector(vector<int>& A)"
A = [1,-2,3,4]
Alst = Py2Cpp_int(A)
show_vector(Alst)

print "\nTest3: int is_in_vector(int a, vector<int>& A)"
print is_in_vector(3, Alst)

print "\nTest3: int is_in_vector(int a, vector<int>& A, vector<int>& pos)"
pos = Py2Cpp_int([0])
res = is_in_vector(-2, Alst, pos)
print "res = ", res, "pos = ", show_vector(pos)

print "\nTest3a: int is_in_vector(int a, vector<int>& A, vector<int>& pos)"
Alst = Py2Cpp_int([1,1,-2,3,-2,-4,5,6,-2])
pos = Py2Cpp_int([0])
res = is_in_vector(-2, Alst, pos)
print "res = ", res, "pos = ", show_vector(pos)


print "\nTest4: is_repeating(vector<int>& A,int& reap)"
Alst = Py2Cpp_int([1,-1,-2,3,-2,-4,5,6,-2])
pos = Py2Cpp_int([0])
print is_repeating(Alst)


print "\nTest5: "
Alst = Py2Cpp_int([1,-1,2,-2])
Blst = Py2Cpp_int([1,-1,2,-3])
print delta(Alst,Blst)


print "\nTest5b: "
Alst = Py2Cpp_int([1,-1,2,-2])
Blst = Py2Cpp_int([1,-1,3,-3])
print delta(Alst,Blst)

print "\nTest5c: "
Alst = Py2Cpp_int([1,-3,2,-2])
Blst = Py2Cpp_int([1,-3,4,-2])
print delta(Alst,Blst)








