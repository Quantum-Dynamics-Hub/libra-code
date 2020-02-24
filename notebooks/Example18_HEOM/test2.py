
import os
import sys
import math
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import * #libdyn
from libra_py import units

nquant = 2
KK = 0

d= nquant * (KK+1) #3
max_tier=2
verbosity=2


vec = intList2()
vec_plus = intList2()
vec_minus = intList2()

gen_hierarchy(d, max_tier, verbosity, vec, vec_plus, vec_minus)



#=============== Previous conventions ==============

LL = max_tier

nn_tot = compute_nn_tot(nquant, KK, LL)
print(F"nn_tot = {nn_tot}")

nn = allocate_3D(nquant+1, KK+1, nn_tot+1)
map_nplus = allocate_3D(nquant+1, KK+1, nn_tot+1)
map_nneg = allocate_3D(nquant+1, KK+1, nn_tot+1)
zero = allocate_1D(nn_tot+1)
map_sum = allocate_1D(LL+1)

compute_nn(nquant, KK, LL, map_sum, nn);

print("nn = \n")
for n in range(1,nn_tot+1):
    print(F"============ ADM index = {n-1} =================")

    line = F"["
    for m in range(1, nquant+1):
        for k in range(KK+1):
            line = line + F" {nn[m][k][n]}," 
    line = line + "]"
    print(line)

compute_map(nquant, KK, LL, nn, map_nplus, map_nneg);

print("map_nplus = \n")
for n in range(1, nn_tot+1):
    print(F"============ ADM index = {n-1} =================")

    line = F" ["
    for m in range(1,nquant+1):
        for k in range(KK+1):
            line = line + F" {map_nplus[m][k][n] - 1}, "
    line = line + "]"
    print(line)

print("map_nneg = \n")
for n in range(1, nn_tot+1):
    print(F"============ ADM index = {n-1} =================")

    line = F" ["
    for m in range(1, nquant+1):
        for k in range(KK+1):
            line = line + F" {map_nneg[m][k][n] - 1}, "
    line = line + "]"
    print(line)




