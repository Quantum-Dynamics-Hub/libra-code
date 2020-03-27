import os
import sys
import math
import time

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import * #libdyn
from libra_py import units


nquant = 2
KK = 1
d= nquant * (KK+1)
max_tier = 12


nn_tot = compute_nn_tot(d, max_tier)
print(F"The expected number of ADMS = {nn_tot} ") 


verbosity=2

vec = intList2()
vec_plus = intList2()
vec_minus = intList2()

t1 = time.time()
gen_hierarchy(d, max_tier, verbosity, vec, vec_plus, vec_minus)
t2 = time.time()

print(F"Time = {t2-t1}")
print(F"Actual number of ADMs = {len(vec)} ") 

