#*********************************************************************************
#* Copyright (C) 2018 Brendan Smith and Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import cmath
import math
import os
import sys
import unittest

cwd = os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/../../_build/src/opt")
sys.path.insert(1,cwd+"/../../_build/src/converters")
sys.path.insert(1,cwd+"/../../_build/src/math_linalg")

# Fisrt, we add the location of the library to test to the PYTHON path
from liblibra_core import *
from libconverters import *
#from libopt import *
#from liblinalg import *

def error_function(q, params):
    """
    q - is a MATRIX(Ndof, 3)
    params - defines the DOFs which we are optimizing

    Now, we just need to define this function
    """

    Ndof = q.num_of_rows     
    active_dofs = params["active_dofs"]  # set grad to 0.0 for any DOFs other than active ones
    k = params["force_constant"]
    d = params["d"]

    grd = MATRIX(Ndof,3) 
  
    # Now, define the gradient for the active dofs

    grad_x, grad_y, grad_z = 0.0, 0.0, 0.0    


    # The repulsion terms between all new atoms and between the new atoms
    # and the old ones bound to the central atoms.
    if params["opt_type"] == 1:
        for i in active_dofs:
            for j in xrange(1,Ndof):
                if i!=j:

                    xi, yi, zi = q.get(i,0) , q.get(i,1) , q.get(i,2)
                    xj, yj, zj = q.get(j,0) , q.get(j,1) , q.get(j,2)
                    a_ij = ( (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2 )**2

                    grad_x += 4.0 * (xi-xj)/a_ij
                    grad_y += 4.0 * (yi-yj)/a_ij
                    grad_z += 4.0 * (zi-zj)/a_ij


    if params["opt_type"] == 2:
        for i in active_dofs:

            xi, yi, zi = q.get(i,0) , q.get(i,1) , q.get(i,2)

            # Harmonic potential terms for the interactions between the new atoms
            # and the central one (this is the atom with index 0)
            x0, y0, z0 = q.get(0,0) , q.get(0,1) , q.get(0,2)
            a_i0 = ( (xi - x0)**2 + (yi - y0)**2 + (zi - z0)**2 )

            grad_x += k*(a_i0 - d)*(xi-x0)/a_i0
            grad_y += k*(a_i0 - d)*(yi-y0)/a_i0
            grad_z += k*(a_i0 - d)*(zi-z0)/a_i0
 
        grd.set(i,0,grad_x)                                 
        grd.set(i,1,grad_y)
        grd.set(i,2,grad_z)

        
    return grd



def convert2matrix(R_central, R_bound, R_new):

    N_bound = len(R_bound) 
    N_new = len(R_new)

  
    X = MATRIX(1+N_bound+N_new, 3)

    X.set(0, 0, R_central.x);   X.set(0, 1, R_central.y);  X.set(0, 2, R_central.z)

    for i in xrange(N_bound):
        X.set(1+i, 0, R_bound[i].x);   X.set(1+i, 1, R_bound[i].y);  X.set(1+i, 2, R_bound[i].z)

    for i in xrange(N_new):
        X.set(1+N_bound+i, 0, R_new[i].x);   X.set(1+N_bound+i, 1, R_new[i].y);  X.set(1+N_bound+i, 2, R_new[i].z)

    return X


def extract_from_matrix(q,N_new):
    """
    Extract the last N_new rows from the matrix and represent them as a list of VECTORS
    """

    R = []
    for i in range(q.num_of_rows-N_new, q.num_of_rows):
        r = VECTOR(q.get(i,0), q.get(i,1), q.get(i,2))
        R.append(r)

    return R

    

def spherical2Cart(r, theta, phi):
    """
    R - coordinates of capping atom in polar coordinates
    """

    R_cart = VECTOR()

    R_cart.x = r*math.sin(theta)*math.cos(phi)
    R_cart.y = r*math.sin(theta)*math.sin(phi)
    R_cart.z = r*math.cos(theta)
   
    return R_cart



        

def cap_atom(R_cent, R_bound, ncap, rnd):
    """
    R_cent - the coordinates of the atom that needs capping (central atom)
    R_bound - the coordinates of the atoms that are bound to the central one
    ncap - the number of the capping atoms needed by the central atom  
    rnd - random numbers generator class

 
    Return: the new set of coordinates R_new - for instance 2 H atoms
    """
    
    # Determine the minial distance between the cetral atom (to be capped)
    # and the atoms that are bound to it    
    bond_len = []
    x_len, y_len, z_len = 0.0, 0.0, 0.0
 
    for i in xrange(len(R_bound)):
        x = R_bound[i].x - R_cent.x #- R_bound[i].x 
        y = R_bound[i].y - R_cent.y #- R_bound[i].y 
        z = R_bound[i].z - R_cent.z #- R_bound[i].z 

        bond_len.append( math.sqrt( x*x + y*y + z*z ) )
    r = 0.5*min(bond_len) 
    d = sum(bond_len)/len(bond_len)

    # Now, add the atoms around the one that needs capping
    R_new = []
    for i in xrange(ncap):
        theta = rnd.uniform(0.0,2.0*math.pi)
        phi = rnd.uniform(0.0,2.0*math.pi)                

        R_new.append( spherical2Cart( r, theta, phi ) + R_cent )   


    # Convert R[n], R_bound, and R_new to the MATRIX type:
    q = convert2matrix(R_cent, R_bound, R_new)

    # Optimize the positions of the added atoms
    opt_tol, opt_step_size, opt_max_steps = 1e-8, 0.01, 1000

    params = {"active_dofs":range(1+len(R_bound), 1+len(R_bound)+len(R_new) ) , "opt_type":1, "force_constant": 1.0, "d": d } 
    q_opt = grad_descent(error_function, q, params, opt_tol, opt_step_size, opt_max_steps)

    params = {"active_dofs":range(1+len(R_bound), 1+len(R_bound)+len(R_new) ) , "opt_type":2, "force_constant": 1.0, "d": d }
    q_opt = grad_descent(error_function, q, params, opt_tol, opt_step_size, opt_max_steps)

    #sys.exit(0)

    R_new = extract_from_matrix(q_opt,len(R_new))  # get the optimized coordinates of the added atoms
    
    return R_new


def cap_system(L, R, connectivity, MaxCoord, new_L):
    """
    L - List of labels for all atoms
    R - List of the coordinates for all atoms
    connectivity - the list of 2-element lists [i, [indiced of the atoms connected to i]]     
    MaxCoord - maximal coordination number for all atoms
    new_L - maps a string input (label of atom to be capped), with a string output (label of capping atom)

    Return - List of labels and coordinates for capped system

    """  

    rnd = Random()
  
    Natoms = len(L)

    L_final = [] 
    R_final = []
    # Copy the original atoms to the results
    for n in xrange(Natoms):
        L_final.append(L[n])
        R_final.append(R[n])

    # For all atoms
    for n in xrange(Natoms):           

        # Find the connectivities, compare with the max coordination
        # and compute the number of valencies to saturate
        neighbors = connectivity[n][1]
        ncaps = MaxCoord[n] - len(neighbors) 

        # Need capping
        if ncaps>0:

            # Collect the information on the local geometry around the
            # atom to be capped
            R_bound = []   # coordinates of all atoms connected to the central one
            for i in neighbors: 
                R_bound.append( R[i] )   

            # Add atoms nearby at random
            R_new = cap_atom(R[n], R_bound, ncaps, rnd)

            # Determine which atom type to add to the "central" one 
            # and add these atoms
            L_new = new_L[ L[n] ]
            for r in R_new:
                L_final.append(L_new)
                R_final.append(r)

                          
    return L_final, R_final


#################################################################################
