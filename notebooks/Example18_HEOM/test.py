
import os
import sys
import math
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import * #libdyn
from libra_py import units



def gen_next_level(parent):
    """
    parent = ( n10, n11, ... , n1K    ...     nM0, nM1, ... , nMK)

    """

    d = len(parent)  # shall be M * (K+1)

    children = []
    for i in range(d):
        children.append( list(parent) )
        children[i][i] += 1

    return children


def gen_next_level2(parents):

    next_level = []

    for parent in parents:
    
        children = gen_next_level(parent)

        for child in children:
            next_level.append(child)

    return next_level



def gen_hierarchy(d, max_tier):
    """
    #    d = nquant * (KK+1)
    # max_tier  = max_level+1
    """
 

    all_vectors = []
    all_coordinates = []  # for each vector - the coordinates are (L, i)
    tier_nums = [] # the number of the nodes for the tier up to a given one


    parents = [ [0]*d ]

    for tier in range(max_tier+1):

        iparent = 0
        for parent in parents:
            if parent not in all_vectors:
                all_vectors.append( list(parent) )
                all_coordinates.append( [tier, iparent] )
                iparent += 1

        tier_nums.append(len(all_vectors)) 

        new_parents = gen_next_level2(parents)
        parents = list(new_parents)


    #===============================================

    print(tier_nums)

    num_nodes = len(all_vectors)


    vec_plus  = []
    vec_minus = []

    for i in range(num_nodes):
        tmp = []
        for j in range(d):
            tmp.append(-1)
        vec_plus.append( list(tmp) )
        vec_minus.append( list(tmp) )


    for i in range(num_nodes):
        L, ipos = all_coordinates[i][0], all_coordinates[i][1]

        for k in range(d):
            np = list(all_vectors[i]); np[k] += 1 

            max_range = max_tier
            if L<max_tier:
                max_range = tier_nums[L+1]

            for j in range(tier_nums[L], max_range):
                if np == all_vectors[j]:
                    vec_plus[i][k] = j


        for k in range(d):
            nm = list(all_vectors[i]); nm[k] -= 1

            min_range = 0
            if L>=2:
                min_range = tier_nums[L-2]
            for j in range(min_range, tier_nums[L-1]):
                if nm == all_vectors[j]:
                    vec_minus[i][k] = j



    for i in range(num_nodes):
        print(F"i= {i}  vec_plus={vec_plus[i]}   vec_minus={vec_minus[i]}")


    return all_vectors, all_coordinates, vec_plus, vec_minus


def print_hierarchy(res, coordinates, vec_plus, vec_minus):

    print(F"size = {len(res)}")
    sz = len(res)

    for i in range(sz):
        print(F"index={i}  vector={res[i]}  coordinate={coordinates[i]} vec_plus={vec_plus[i]} vec_minus={vec_minus[i]}")

    
res, coords, vec_plus, vec_minus = gen_hierarchy(3, 3)
print_hierarchy(res, coords, vec_plus, vec_minus)




nquant = 2
LL = 10
KK = 0


nn_tot = compute_nn_tot(nquant, KK, LL)
print(F"nn_tot = {nn_tot}")

res, coords, vec_plus, vec_minus = gen_hierarchy(nquant*(KK+1), LL)
print_hierarchy(res, coords, vec_plus, vec_minus)


sys.exit(0)


nn = allocate_3D(nquant+1, KK+1, nn_tot+1)
map_nplus = allocate_3D(nquant+1, KK+1, nn_tot+1)
map_nneg = allocate_3D(nquant+1, KK+1, nn_tot+1)
zero = allocate_1D(nn_tot+1)
map_sum = allocate_1D(LL+1)

compute_nn(nquant, KK, LL, map_sum, nn);

print("nn = \n")
for L in range(nn_tot+1):
    print(F"============ L = {L} =================")
    for n in range(nquant+1):
        line = F"n = {n} ["
        for k in range(KK+1):
            line = line + F" {nn[n][k][L]}"
        line = line + "]"
        print(line)

compute_map(nquant, KK, LL, nn, map_nplus, map_nneg);

print("map_nplus = \n")
for L in range(nn_tot+1):
    print(F"============ L = {L} =================")
    for n in range(nquant+1):
        line = F"n = {n} ["
        for k in range(KK+1):
            line = line + F" {map_nplus[n][k][L]}"
        line = line + "]"
        print(line)

print("map_nneg = \n")
for L in range(nn_tot+1):
    print(F"============ L = {L} =================")
    for n in range(nquant+1):
        line = F"n = {n} ["
        for k in range(KK+1):
            line = line + F" {map_nneg[n][k][L]}"
        line = line + "]"
        print(line)


