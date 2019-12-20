import os
import sys
import math
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import * #libdyn
from libra_py import units



nquant = 2
LL = 10
KK = 0


nn_tot = compute_nn_tot(nquant, KK, LL)
print(F"nn_tot = {nn_tot}")


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




def unpack_rho_py(RHO):

    res = CMATRIXList()
    rho = CMATRIX(nquant, nquant)

    sub_y = Py2Cpp_int(list(range(nquant))) 
    
    for n in range(0, nn_tot+1):

        sub_x = Py2Cpp_int(list(range(n*nquant, (n+1)*nquant)))        
        pop_submatrix(RHO, rho, sub_x, sub_y )

        res.append( CMATRIX(rho) )

    return res


def pack_rho_py(rho_unpacked, RHO):

    sub_y = Py2Cpp_int(list(range(nquant))) 
    
    for n in range(0, nn_tot+1):
        sub_x = Py2Cpp_int(list(range(n*nquant, (n+1)*nquant)))        
        push_submatrix(RHO, rho_unpacked[n], sub_x, sub_y )


    
    

def compute_derivatives(RHO, prms):
    """
    Here we simply return the pre-computed global variables
    """
    Ham = prms["Ham"]    
    
    dRHO = CMATRIX((nn_tot+1)*nquant, nquant)  # all drho_dt matrices stacked on top of each other
    rho = CMATRIX(nquant, nquant)
    drho = CMATRIX(nquant, nquant)

    sub_y = Py2Cpp_int(list(range(nquant))) 


    rho_unpacked = unpack_rho_py(RHO)

    for n in range(1, nn_tot+1):

        sub_x = Py2Cpp_int(list(range(n*nquant, (n+1)*nquant)))
        
        drho = compute_deriv_n(n, rho_unpacked, prms["Ham"], prms["eta"], prms["temperature"], prms["gamma_matsubara"], prms["c_matsubara"],
                               prms["nn"], prms["KK"], prms["zero"], prms["map_nplus"],  prms["map_nneg"])

        push_submatrix(dRHO, drho, sub_x, sub_y )

    return dRHO






params = {}

#============== HEOM topology ==============


params.update( {  "KK":KK,
                  "nn":nn, 
                  "zero":zero, 
                  "map_nplus":map_nplus,
                  "map_nneg":map_nneg 
                } )


#============== System ==============
Ham = CMATRIX(2,2)
Ham.set(0, 0, 50.0 * units.inv_cm2Ha);    Ham.set(0, 1, 200.0 * units.inv_cm2Ha);
Ham.set(1, 0,200.0 * units.inv_cm2Ha);    Ham.set(1, 1, -50.0 * units.inv_cm2Ha);

params.update( { "Ham": Ham } )


#============== Bath ==============
params.update({ "gamma": 1.0/(0.1 * units.ps2au),
                "eta": 2.0 * 50.0 * units.inv_cm2Ha,
                "temperature": 300.0 
             })

gamma_matsubara = doubleList()
c_matsubara = complexList()
setup_bath(params, gamma_matsubara, c_matsubara)

for k in range(KK+1):
    print(F" k = {k} gamma_matsubara[{k}] = {gamma_matsubara[k]}  c_matsubara[{k}] = {c_matsubara[k]}")

params.update({ "gamma_matsubara": gamma_matsubara, "c_matsubara":c_matsubara  } )



#============= Initialization ============

rho_unpacked = CMATRIXList()
for n in range(nn_tot+1):
    rho_unpacked.append( CMATRIX(nquant, nquant))

rho = CMATRIX((nn_tot+1)*nquant, nquant)  # all rho matrices stacked on top of each other
drho = CMATRIX((nn_tot+1)*nquant, nquant)  # all drho_dt matrices stacked on top of each other


# Initial conditions
rho.set(2, 0, 1.0+0.0j)  # rho[1](0,0)


#============== Propagation =============


f = open("pops.txt", "w")
f.close()

dt = 0.1 * units.fs2au
nsteps = 10000
tolerance = 1e-6


for step in range(nsteps):
    time = dt * step / units.fs2au

    f = open("pops.txt", "a")
    f.write( " %5.3f  %10.8f  %10.8f\n" % (time, rho.get(2,0).real, rho.get(3,1).real))
    f.close()

    if step % 10 == 1:
        unpack_rho(rho_unpacked, rho)
        params["zero"] = filter(rho_unpacked, tolerance);
        pack_rho(rho_unpacked, rho)

    rho = RK4(rho, dt, compute_heom_derivatives, params)






