#*********************************************************************************
#* Copyright (C) 2021 Mohammad Shakiba, Alexey V. Akimov
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 3 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
"""
.. module:: molden_methods
   :platform: Unix, Windows
   :synopsis: This module implements functions for processing the molden files printed out for molecular orbitals (works properly for CP2K outputs).
.. moduleauthors:: 
       Mohammad Shakiba, Alexey V. Akimov 
  
"""

import os
import sys
import numpy as np
import time

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import data_outs
from libra_py import units
import math


def molden_file_to_libint_shell(molden_filename: str, is_spherical: bool, is_periodic=False,
                                cell=np.array([[0,0,0],[0,0,0],[0,0,0]]), R_vec=np.array([0,0,0])):
    """
    This function gets the molden file and returns the shell for use with 
    libint to compute the atomic orbital overlaps.
    
    Args:
    
        molden_filename (string): The name of the molden file.
        
        is_spherical (bool): The Gaussian cartesian or spherical orbital flag.
        
    Returns:
    
        shell (shell data type from libint): The shell containing the exponents and contraction  
                                             coefficients for each atom, their coordinates in Bohr,
                                             and their angular momentum values.
        l_vals_all (list): The angular momentum values for all the atoms and basis set in order.
                           This is used later to resort the eigenvectors in molden file.
        
    """
    # Reading the lines of the molden file
    file = open(molden_filename,'r')
    lines = file.readlines()
    file.close()

    # The characters in the molden file start with '[' including
    # coordinates, basis sets, eigenvectors, etc. We store them
    # into chars to move between them in a for loop
    chars = []
    for i in range(len(lines)):
        if lines[i].replace(' ','')[0]=='[':
            chars.append(i)
    
    # The coordinates of the structure, these are in Bohr value
    coords = []
    for i in range(chars[1]+1,chars[2]):
        # Each line in the coordinates part
        tmp = lines[i].split()
        # not necessary but let's keep it, might be useful in future
        atomic_number = int(tmp[2])
        # The x, y, and z coordinates
        x = float(tmp[3])
        y = float(tmp[4])
        z = float(tmp[5])
        if is_periodic:
            translation = np.dot(R_vec,cell)
            x += translation[0]
            y += translation[1]
            z += translation[2]
        # append all in the coords
        coords.append([atomic_number, x, y, z])
        
    
    """
    The molden format for basis set is the same as below:
     [GTO]
        1       0        * this line
                         s       6    1.00 * three_line
                                                       106.389717      -0.004151
                                                        19.510792      -0.020670
                                                         5.482954      -0.051503
                                                         0.784075       0.334627
                                                         0.355861       0.562106
                                                         0.169708       0.171299
                         p       6    1.00 * three_line
                                                        19.697346       0.007924
                                                         5.136667       0.051441
                                                         1.837949       0.189840
                                                         0.768298       0.404986
                                                         0.351318       0.401236
                                                         0.166091       0.105186
    
        2       0        * atom_line
                         s       6    1.00 * three_line
                                                       106.389717      -0.004151
                                                        19.510792      -0.020670
                                                         5.482954      -0.051503
                                                         0.784075       0.334627
                                                         0.355861       0.562106
                                                         0.169708       0.171299
                         p       6    1.00 * three_line
                                                        19.697346       0.007924
                                                         5.136667       0.051441
                                                         1.837949       0.189840
                                                         0.768298       0.404986
                                                         0.351318       0.401236
                                                         0.166091       0.105186
    
    What we need to do is to first find the initial lines I showed with a star* which show
    the basis set for each atom number.
    """
    # atom lines
    atom_line = []
    # first need to set up a counter to check whether 
    # it can be an integer or not. We move between [GTO] and [MO]
    # in the file.
    c = 0
    for i in range(chars[2]+1,chars[3]):
        try:
            # the first value, if it can be turned into 
            # an integer value
            first_val = int(lines[i].split()[0])
            # counter starts to increase now.
            c += 1
            # if it was equal to counter add it to the list
            if first_val==c:
                atom_line.append(i)
        except:
            pass
        
    # Setting up the counters
    # the counter for counting the atoms
    # necessary to initialize the shell for the first one
    # then for others we add to the initialized shell 
    # using add_to_shell
    c = 0
    # coords counter
    h1 = 0
    # all the angular momentum values for all atoms
    l_vals_all = []
    for i in range(len(atom_line)):
        tmp = []
        # Lines with three elements in them as is shown above
        three_line = []
        # the angular momentum values for each atom
        l_vals = []
        if i==len(atom_line)-1:
            # next atom_line
            end_l = chars[3]
        else:
            # next atom_line
            end_l = atom_line[i+1]
        # now loop over them and append the angular momentum values
        # for each atom.
        for j in range(atom_line[i],end_l):
            if len(lines[j].split())==3:
                three_line.append(j)
                # if s orbital
                if lines[j].split()[0]=='s':
                    l_vals.append(0)
                    l_vals_all.append(0)
                # if p orbital
                elif lines[j].split()[0]=='p':
                    l_vals.append(1)
                    l_vals_all.append(1)
                # if d orbital
                elif lines[j].split()[0]=='d':
                    l_vals.append(2)
                    l_vals_all.append(2)
                # if f orbital
                elif lines[j].split()[0]=='f':
                    l_vals.append(3)
                    l_vals_all.append(3)
                # if g orbital
                elif lines[j].split()[0]=='g':
                    l_vals.append(4)
                    l_vals_all.append(4)
        # a counter for the l_vals
        c1 = 0
        for j in three_line:
            x = []
            # the number of exponents and contraction coefficients 
            n_line = int(lines[j].split()[1])
            c += 1
            # exponents
            exp = []
            # contraction coefficients
            coeff = []
            # now the loop for the exponents and coefficients
            for k in range(1,n_line+1):
                # the line index
                n = j+k
                # append the exponents and coefficients
                exp.append(float(lines[n].split()[0]))
                coeff.append(float(lines[n].split()[1]))
            # the coordinates of the atom
            coords_init = [coords[h1][1],coords[h1][2],coords[h1][3]]
            # This is how we initialize an Atom data type for the coordinats to use in libint
            a = VECTOR(coords_init[0], coords_init[1], coords_init[2])
            # If the counter is the first one
            if c==1:
                # initialize the shell by creating a shell type object
                shell = initialize_shell(int(l_vals[c1]), is_spherical, Py2Cpp_double(exp), Py2Cpp_double(coeff), a)
            else:
                # we now use the add_to_shell function to add the necessary values to the shell
                add_to_shell(shell, int(l_vals[c1]), is_spherical, Py2Cpp_double(exp), Py2Cpp_double(coeff), a)
            c1 += 1
        h1 += 1
    
    return shell, l_vals_all


def eigenvectors_molden(molden_filename: str, nbasis: int, l_vals: list):
    """
    This functions read the molecular orbitals eigenvectors and resort 
    them so that we can use it with the libint.
    
    Args:
    
        molden_filename (string): The name of molden file.
        
        nbasis (integer): This is obtained from the liblibra_core.nbasis. It is necessary 
                          since sometime the molden files will not append all the values
                          and there would be some values missing. It is better to use a 
                          higher value of NDIGITS in the CP2K input file to have most of
                          the values.
        
        l_vals (list): All the angular momentum values return by the molden_file_to_libint_shell
                       function for the molden file.
                       
    Returns:
        
        eigenvectors (list): All the sorted eigenvectors for use with the libint computed 
                                    atomic orbital overlaps.
        
        orbital_energies (numpy array): All orbital energies for the eigenvectors.
    
    """
    # again reading the molden files and its lines
    file = open(molden_filename,'r')
    lines = file.readlines()
    file.close()
    
    # energy lines
    orbital_energies = []
    ener_lines = []
    for i in range(len(lines)):
        if 'Ene='.lower() in lines[i].lower():
            ener_lines.append(i)
            # removing the 'Ene=' from the lines
            tmp = lines[i].lower().replace('ene=','')
            # the energy value
            energy = float(tmp.split()[0])
            # appending the energy into orbital_energies
            orbital_energies.append(energy)
    # make it a numpy array
    orbital_energies = np.array(orbital_energies)
    
    # eigenvectors list
    eigenvectors = []
    # loop over the ener_lines
    for i in range(len(ener_lines)):
        # With this if there is no value in
        # the molden file it would automatically be zero
        eigenvector = np.zeros((nbasis))
        if i==len((ener_lines))-1:
            # end line to loop over for an energy value
            endl = len(lines)
        else:
            # end line to loop over for an energy value
            endl = ener_lines[i+1]
        # now loop between them
        for j in range(ener_lines[i]+3,endl):
            # the index in the molden file
            c = int(lines[j].split()[0]) - 1
            # append the eigenvector value to 
            # the np.zero eigenvector
            eigenvector[c] = float(lines[j].split()[1])
        
#         eigenvector = np.array(eigenvector)
        # new indices form resot eigenvectors based on l_vals
        new_indices = resort_eigenvectors(l_vals)
        # the new and sorted eigenvector
        eigenvector = eigenvector[new_indices]
        # append it to the eigenvectors list
        eigenvectors.append(eigenvector)
        
    return eigenvectors, orbital_energies

def resort_eigenvectors(l_vals):
    """
    This function resorts the indices of the eigenvectors based on the 
    angular momentum values. 
    
    Args:
    
        l_vals (list): All the values of the angular momentum for all the atoms.
        
    Returns:
    
        new_indices (list): The new indices.
        
    """
    # new indices
    new_indices = []
    
    """
    The difference between the way MOLog and molden files 
    store the eigenvectors is like this (as an example for the Cd atom):
    
    MOLog files        Molden files
 
    Cd  2s            Cd  2s   
    Cd  3s            Cd  3s   
    Cd  3py           Cd  3px
    Cd  3pz           Cd  3py  
    Cd  3px           Cd  3pz  
    Cd  4py           Cd  4px  
    Cd  4pz           Cd  4py  
    Cd  4px           Cd  4pz  
    Cd  4d-2          Cd  4d0
    Cd  4d-1          Cd  4d+1
    Cd  4d0           Cd  4d-1 
    Cd  4d+1          Cd  4d+2 
    Cd  4d+2          Cd  4d-2 
    Cd  5d-2          Cd  5d0 
    Cd  5d-1          Cd  5d+1 
    Cd  5d0           Cd  5d-1 
    Cd  5d+1          Cd  5d+2 
    Cd  5d+2          Cd  5d-2 
    Cd  5f-3          Cd  5f0 
    Cd  5f-2          Cd  5f+1
    Cd  5f-1          Cd  5f-1 
    Cd  5f0           Cd  5f+2 
    Cd  5f+1          Cd  5f-2 
    Cd  5f+2          Cd  5f+3 
    Cd  5f+3          Cd  5f-3 
    """
    
    # setting up the counter
    c = 0
    # loop over all the angular momentum values
    for i in range(len(l_vals)):
        l_val = l_vals[i]
        # find the reorder indices needed for this l_val
        reordered_indices = index_reorder(l_val)
        for j in range(len(reordered_indices)):
            # now append it by plus the counter since
            # we aim to do it for all the eigenvectors
            # and l values
            new_indices.append(c+reordered_indices[j])
        # increase the counter with t
        c += len(reordered_indices)
        
    # Return the new indices
    return new_indices
    
    
def index_reorder(l_val):
    """
    This function returns the reordered indices for an angular momentum value.
    
    Args:
    
        l_val (integer): The angular momentum value.
        
    Returns:
        
        new_order (numpy array): Contains the new indices.
    
    """
    # for s orbital
    if l_val == 0:
        new_order = [1]
    # for p orbital
    elif l_val == 1:
        new_order = [2,3,1]
    # for d orbital
    elif l_val == 2:
        new_order = [5,3,1,2,4]
    # for f orbital
    elif l_val == 3:
        new_order = [7,5,3,1,2,4,6]
    # for g orbital
    elif l_val == 4:
        new_order = [9,7,5,3,1,2,4,6,8]

    # The indeices
    return np.array(new_order)-1

