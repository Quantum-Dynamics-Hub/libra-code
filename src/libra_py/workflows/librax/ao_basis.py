#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file ao_basis.py
# This module implements the functions that constructs atomic orbital basis
# by n Gaussian Type Orbitals (nGTO)

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def input_AO_name(l_atoms, atom_spec, basis_type, flag):
    ##
    # \param[in] l_atoms A list of all atom labels (e.g. C, C, C, H, H, ..)
    # \param[in] atom_spec A list of distinct atomic species (e.g. C, H)
    # \param[in] basis_type A list of AO types coming with the atom of each distinct type
    # (e.g. basis_type[i] = ["S","P"] implies that the atom of type i (e.g. C) has S and P orbitals) 
    # \param[in] flag Debug info printing: 1 - print, otherwise - don't
    #
    # This function returns the list of atomic orbital type (s, px, py, pz, etc...) and the number
    # of orbitals of each type: S, P, D, L (S + P)
    #
    # Used in: extract_gms(g09).py/gms(g09)_extract_ao_basis

    orb_name = []

    for la in l_atoms: # all atoms

        # Determine the kind of the atom
        i = -1
        for j in range(0,len(atom_spec)): 
            if la == atom_spec[j]:
                i = j
        if i==-1:
            print "Error in input_AO_name: atom ", la, " is not known\n"
            sys.exit(0)

        aoa_tmp = []
        orb_name1 = []
        scount = 0
        pcount = 0
        dcount = 0
        lcount = 0

        # Over all basis functions coming with the atom of type i
        # len(basis_type[i]) - the number of AOs coming with atom of type i
        for j in range(0,len(basis_type[i])): 
            b_tmp = basis_type[i][j]
            if b_tmp == "S":
                scount += 1
                if scount < 2:
                    orb_name1.append("s")
            if b_tmp == "P":
                pcount += 1
                if pcount < 2:
                    orb_name1.append("px")
                    orb_name1.append("py")
                    orb_name1.append("pz")
            if b_tmp == "D":
                dcount += 1
                if dcount < 2:
                    orb_name1.append("dxy")
                    orb_name1.append("dyz")
                    orb_name1.append("dzx")
                    orb_name1.append("dx^2-y^2")
                    orb_name1.append("dz^2")
            elif b_tmp == "L":
                lcount += 1
                if lcount < 2:
                    orb_name1.append("s")
                    orb_name1.append("px")
                    orb_name1.append("py")
                    orb_name1.append("pz")
        orb_name.append(orb_name1)


    if flag == 1:
        print "number of s-type orbitals =", scount
        print "number of p-type orbitals =", pcount
        print "number of d-type orbitals =", dcount
        print "number of l-type orbitals =", lcount
        print "all orbitals: "
        for i in xrange(len(orb_name)):
            print i, orb_name[i]


    return orb_name, scount, pcount, dcount, lcount



def construct_ao_basis(ao_data,label,R,nGTO,orb_name):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : The list which contains extracted data from the file.
    # This function returns the list of atomic orbital basis as "ao_basis".
    #
    # Used in: extract.py/extract_ao_basis

    expo_s = ao_data["expo_s"]
    expo_p = ao_data["expo_p"]
    expo_d = ao_data["expo_d"]
    coef_s = ao_data["coef_s"]
    coef_p = ao_data["coef_p"]
    coef_d = ao_data["coef_d"]

    ao_basis = []

    Bohr_to_Angs = 0.529177208
    # define atomic orbitals by using the Gaussian basis
    for i in range(0,len(label)): # all atoms

        k_s = 0 # add nGTO for s orbital
        k_p = 0 #              p orbital

        for j in range(0,len(orb_name[i])):

            ao = AO()  # this is the AO we create - the j-th AO on the i-th atom

            if orb_name[i][j][0] == "s":
                expo_tmp = expo_s[i][k_s:k_s+nGTO];   coef_tmp = coef_s[i][k_s:k_s+nGTO]
                nx,ny,nz = 0,0,0
                k_s += nGTO
            elif orb_name[i][j][0] == "p":
                expo_tmp = expo_p[i][k_p:k_p+nGTO];   coef_tmp = coef_p[i][k_p:k_p+nGTO]
                if orb_name[i][j][1] == "x":
                    nx,ny,nz = 1,0,0
                elif orb_name[i][j][1] == "y":
                    nx,ny,nz = 0,1,0
                elif orb_name[i][j][1] == "z":
                    nx,ny,nz = 0,0,1
                    k_p += nGTO

            elif orb_name[i][j][0] == "d":
                expo_tmp = expo_p[i][k_p:k_p+nGTO];   coef_tmp = coef_p[i][k_p:k_p+nGTO]
                if orb_name[i][j][1:3] == "xy":
                    nx,ny,nz = 1,1,0
                elif orb_name[i][j][1:3] == "yz":
                    nx,ny,nz = 0,1,1
                elif orb_name[i][j][1:3] == "zx":
                    nx,ny,nz = 1,0,1
                elif orb_name[i][j][1:3] == "z^":  # !!!! Here is unfixed BUG !!! 
                    nx,ny,nz = 0,0,2
                # in the case of dx^2-y^2, i should add AO later. 

            # Construct AO from the primitive Gaussians
            for k in range(0,len(expo_tmp)):         
                g = PrimitiveG(nx, ny, nz, expo_tmp[k], R[i]) # single point

                # Contraction coefficients correspond to the Gaussian primitives as they are
                ao.add_primitive(coef_tmp[k], g )
                   
            # Normalize the overall contraction
            ao.normalize()
            ao_basis.append(ao)
    
    return ao_basis


