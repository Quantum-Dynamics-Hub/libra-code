#*********************************************************************************
#* Copyright (C) 2017 Kosuke Sato, Alexey V. Akimov
#* Olga S. Bokareva
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file extract_g09.py
# This module implements the functions that extract
# atomic forces , molecular energies, molecular orbitals, and atomic basis information
# written in Gaussian output file.

import os, re
import sys
import math
import itertools as it
import numpy as np

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import detect_g09
import ao_basis

Q_l_dict = {1.0 : 'H', 6.0 : 'C', 7.0 : 'N', 8.0 : 'O'}

def g09_extract_ao_basis(inp_str, label, R, flag): # DONE!!!
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] inp_str A list containing info. for atomic orbital basis.
    # \param[in] label   A list of atomic labels (e.g. C,H,O)
    # \param[in] R       A list of atomic coordinates (R[i] =[R[i].x, R[i].y, R[i].z]
    # \param[in] flag    debugging info: option 1-> print, otherwise -> don't 
    # ao - returned lists of atomic orbital basis sets 
    #
    # Used in: extract_g09.py/g09_extract

    # atomic species
    l_atom_spec = []
    atom_spec = []
    sz = len(inp_str)

    # atom labels should be taken from the coordinates because they are not given in the first basis set section
    atom_spec = label

    for i in xrange(sz):
        spline = inp_str[i].split()
        if len(spline) == 2 and spline[1] == '0':
#    	    print spline
            l_atom_spec.append(i)

    # atomic basis sets
    basis_type = []
    basis_expo = []
    basis_coef = []
    
#    print "l_atom_spec is ", l_atom_spec

    for i in range(0,len(atom_spec)):
        stmp = atom_spec[i]
        type_tmp = []
        expo_tmp = []
        coef_tmp = []
        if i < len(atom_spec) - 1:        # if not the last atom is considered...
            tmp_start = l_atom_spec[i] + 1
            tmp_end = l_atom_spec[i+1] - 2
        else: 			          # if the last atom is considered
            tmp_start = l_atom_spec[i] + 1
            tmp_end = sz - 1                 #ab_end

# 	Reading basis from G09 format
	functions = []
	for j in range(tmp_start,tmp_end+1):
    	    spline = inp_str[j].split()
	    if len(spline)== 4:
		functions.append(spline[0])
		inp_str[j]="*****"

	transformed=[]
	for j in range(tmp_start,tmp_end+1):
	    transformed.append(inp_str[j].strip())
	
	groups = []
	for key,group in it.groupby(transformed,lambda line: line.startswith('*****')):
	    if not key:
	        group = list(group)
	        groups.append(group)    	
	
	for x in xrange(0,len(functions)):
	    for y in xrange(0, len(groups[x])):
		type_tmp.append(functions[x])
		sg = groups[x][y].split()
		coef_tmp1 = []

		if len(sg) == 2:
        	    var0 = re.sub(r'D', 'e', sg[0])
		    var1 = re.sub(r'D', 'e', sg[1])
        	    expo_tmp.append(float(var0))
        	    coef_tmp1.append(float(var1))
        	if len(sg)==3:                        # should be L type of basis functions
        	    var0 = re.sub(r'D', 'e', sg[0])
        	    var1 = re.sub(r'D', 'e', sg[1])
        	    var2 = re.sub(r'D', 'e', sg[2])
        	    expo_tmp.append(float(var0))
        	    coef_tmp1.append(float(var1))
        	    coef_tmp1.append(float(var2))

    		coef_tmp.append(coef_tmp1)
	
	coef_tmp2 = [x for x in coef_tmp if x != []]
        
        basis_type.append(type_tmp)
        basis_expo.append(expo_tmp)
        basis_coef.append(coef_tmp2)
		

    print "basis_type", basis_type
    print "basis_expo", basis_expo
    print "basis_coef", basis_coef

    # ******* get the Gaussian primitives parameters *******
    expo_s = []
    expo_p = []
    expo_d = []
    coef_s = []
    coef_p = []
    coef_d = []

    for la in label: # all atoms
        for j in range(0,len(atom_spec)): # specify the kind of the atom
            if la == atom_spec[j]:
                i = j

        expo_stmp = []
        expo_ptmp = []
        expo_dtmp = []

        coef_stmp = []
        coef_ptmp = []
        coef_dtmp = []

        for j in range(0,len(basis_type[i])): # basis number of atoms
            b_tmp = basis_type[i][j]
#            print "basis type is ", b_tmp
            if b_tmp == "S":
#        	print basis_expo[i][j]
#        	print basis_coef[i][j]
                expo_stmp.append(basis_expo[i][j])
                coef_stmp.append(basis_coef[i][j][0])
            elif b_tmp == "P":
                expo_ptmp.append(basis_expo[i][j])
                coef_ptmp.append(basis_coef[i][j][0])
            elif b_tmp == "D":
                expo_dtmp.append(basis_expo[i][j])
                coef_dtmp.append(basis_coef[i][j][0])
            elif b_tmp == "SP": 			# In G09, SP is used for notation of L basis functions
                expo_stmp.append(basis_expo[i][j])
                expo_ptmp.append(basis_expo[i][j])
#                print "basis_expo are ", basis_expo[i][j]
#                print "basis coefs are ", basis_coef[i][j]
                coef_stmp.append(basis_coef[i][j][0])
                coef_ptmp.append(basis_coef[i][j][1])
            # f orbitals are not taken into account, so should add them.

        expo_s.append(expo_stmp)
        expo_p.append(expo_ptmp)
        expo_d.append(expo_dtmp)
        coef_s.append(coef_stmp)
        coef_p.append(coef_ptmp)
        coef_d.append(coef_dtmp)


    orb_name, scount, pcount, dcount, lcount = ao_basis.input_AO_name(label, atom_spec, basis_type, flag)

    ao_data = {}
    ao_data["expo_s"] = expo_s
    ao_data["expo_p"] = expo_p
    ao_data["expo_d"] = expo_d
    ao_data["coef_s"] = coef_s
    ao_data["coef_p"] = coef_p
    ao_data["coef_d"] = coef_d

    if flag == 1:
        print "expo_s=",ao_data["expo_s"]
        print "expo_p=",ao_data["expo_p"]
        print "expo_d=",ao_data["expo_d"]
        print "coef_s=",ao_data["coef_s"]
        print "coef_p=",ao_data["coef_p"]
        print "coef_d=",ao_data["coef_d"]


    ao = ao_basis.construct_ao_basis(ao_data,label,R,scount,orb_name)

    return ao
    ### This part should work
    
def g09_extract_mo(inp_str,Ngbf,active_space,flag): # DONE
    ##
    # Extracts MO-LCAO coefficients from the the list of input lines
    # 
    # \param[in] inp_str      Strings containing the info for molecular orbitals
    # \param[in] Ngbf         Number of Gaussian Basis Functions
    # \param[in] active_space A list of molecular orbitals considered during the calculation
    # \param[in] flag         debugging info: option 1-> print, otherwise -> don't
    # E - returned MATRIX object, containing the total excitation energies
    # C - returned MATRIX object, containing the eigenvectors:
    # C.get(a,i) - is the coefficient of AO with index a in the MO with index i
    #
    # Used in: extract_g09.py/g09_extract

    stat_span = Ngbf + 3 # period for beginning of coefficient lines

    sz = len(inp_str) 

    # create objects of MATRIX type, containing eigenvalues and eigenvectors 
    E_full = MATRIX(Ngbf,Ngbf)
    C_full = MATRIX(Ngbf,Ngbf)

    for i in xrange(sz):

        if i % stat_span == 0 :
            ind_of_eig = inp_str[i].split() # split lines for indexes of eigenvalues

            # eigenvalues
            eig_val = inp_str[i+2].split() # split lines for eigenvalues
	    eig_val = eig_val[2:] 	   # delete "Eigenvectors --" from g09 output line
	    print "eig_values are", eig_val
            for j in range(0,len(ind_of_eig)):
                k = int(ind_of_eig[j]) - 1           # python index start from 0
                E_full.set(k,k,float(eig_val[j]))

            # molecular coefficients
            ic = i + 3                              # beginning of coefficient lines
            for j in range(ic,ic+Ngbf):             # loop for AO basis
                eig_vec = inp_str[j].split()        # split lines for eigenvectors
                for k in range(0,len(ind_of_eig)):  
                    ja = int(eig_vec[0]) - 1            # index of AO basis
                    ke = int(ind_of_eig[k]) - 1         # index of eigenvectors
		    if eig_vec[1].isdigit():		# if number of atom is in line at second position
			kvec = k + 4
		    else:
			kvec = k + 2
                    C_full.set(ja,ke,float(eig_vec[kvec]))

#    print "E_full is ", E_full.show_matrix()
#    print "C_full is ", C_full.show_matrix()
    
    #Here, we need to add reduction of E, C matrices and the corresponding trunctaion of the ao list. 
    #This will allow us make further computations faster and in consistent with those adopted in QE 

    # ***********Here, reduce E_full and C_full ***************
    sz = len(active_space)
    if sz==0:
        print "active space is not defined correctly, exit....."
        sys.exit(0)

    E = MATRIX(sz,sz)
    for i in xrange(sz):
        imo = active_space[i]-1
        E.set(i,i,E_full.get(imo,imo))

    #C = CMATRIX(C_full) # input full matrix
    C = CMATRIX(Ngbf,sz)
    for i in xrange(Ngbf):
        for j in xrange(sz):
            jmo = active_space[j]-1
            C.set(i,j,C_full.get(i,jmo),0.0)

    # ************************************************************

    if flag == 1:
        print "*** full matrix ****"
        print "E_full(0,0) is",E_full.get(0,0)
        print "E_full(Ngbf-1,Ngbf-1)",E_full.get(Ngbf-1,Ngbf-1)
        print "C_full(0,0) is",C_full.get(0,0)
        print "C_full(Ngbf-1,0) is",C_full.get(Ngbf-1,0)
        print "C_full(0,Ngbf-1) is",C_full.get(0,Ngbf-1)
        print "C_full(Ngbf-1,Ngbf-1) is",C_full.get(Ngbf-1,Ngbf-1)
        print "*** reduced matrix ****"
        print "active_space=",active_space
        print "E Matrix is"; E.show_matrix()
        print "C Matrix is"; C.show_matrix()

    return E, C
    
def g09_extract_coordinates(inp_str,flag): # DONE!
    ##
    # Extracts atomic labels, nuclear charges, and coordinates of all atoms
    # from the the list of input lines
    # each input line is assumed to have the format:
    # number  Q  ....  x y z
    # the label is defined via Q_l_dict listing all atomic names for nuclear charges
    # \param[in] inp_str  Strings containing the info for all atomic coordinates
    # \param[in] flag     debugging info: option 1-> print, otherwise -> don't
    # label - returned list of atomic labels (strings)
    # Q - returned list of nuclear charges (floats)
    # R - returned list of nuclear coordinates (VECTOR objects)
    #
    # Used in: extract_g09.py/g09_extract


    label, Q, R = [], [], []
    for a in inp_str: 
        spline = a.split() 

        # atomic charges
        Q_tmp = float(spline[1])
        Q.append(float(spline[1]))
        label.append(Q_l_dict[Q_tmp])

        # coordinates of atoms: Last 3 elements of the line
        x = float(spline[len(spline)-3])
        y = float(spline[len(spline)-2])
        z = float(spline[len(spline)-1])
        r = VECTOR(x,y,z)
        R.append(r)

    if flag == 1:
        print "label=", label
        print "charges=", Q
        print "coor_atoms="
        for r in R:
            print R.index(r), r, r.x, r.y, r.z


    return label, Q, R


def g09_extract_gradient(inp_str,flag): # DONE!!!
    ##
    # Extracts atomic gradients on all atoms from the the list of input lines
    # each input line is assumed to have the format:
    #  ....  gx  gy  gz
    #
    # \param[in] inp_str  Strings containing the gradient for all atoms
    # \param[in] flag     debugging info: option 1-> print, otherwise -> don't
    # grad - returned list of VECTOR objects
    #
    # Used in: extract_g09.py/g09_extract

    grad = []
    for a in inp_str:
        spline = a.split()  

        # Last 3 elements of the line
        x = float(spline[len(spline)-3])
        y = float(spline[len(spline)-2])
        z = float(spline[len(spline)-1])
        g = VECTOR(-x,-y,-z)  # Note: Gaussian outputs forces
                              # not gradients

        grad.append(g)

    if flag == 1:
        print "atomic gradient="
        for g in grad:
            print grad.index(g), g, g.x, g.y, g.z

    return grad


def g09_extract_first(filename,flag): 
    ##
    # This function only extracts number of electrons in a system.
    # In the main.py, this information is needed to construct active space
    #
    # \param[in] filename    The name of the Gaussian09 output file from which we will be getting info
    # info["Nele"] returns total number of electrons
    # 
    # Used in: main.py/main

    f = open(filename,"r")
    A = f.readlines()
    f.close()

    # detect the lines including information from gamess output file
    info = detect_g09.detect(A,flag)

    return info["Nele"]



def g09_extract(filename,states,min_shift,active_space,flag): # ONLY SLIGHTLY MODIFIED
    ##
    # This function extracts all the necessary information (energies, gradients, coordinates, MOs, 
    # AOs, etc. ) from the Gaussian09 output.
    #
    # \param[in] filename    The name of the Gaussian09 output file from which we will be getting info
    # \param[in] states      excitation states
    # \param[in] min_shift   e.g. -1 -> includes HOMO-1, HOMO
    # \param[in] active_space molecular orbital considered during the calculation 
    # \param[in] flag  debugging info: option 1-> print, otherwise -> don't
    # label - returned list of atomic labels (strings)
    # Q - returned list of nuclear charges (floats)
    # R - returned list of nuclear coordinates (VECTOR objects)
    # grad - returned list of VECTOR objects
    # E - returned MATRIX object, containing the total excitation energies 
    # C - returned MATRIX object, containing the eigenvectors: 
    # ao - returned lists of atomic orbital basis sets 
    #
    # Used in: main.py/main or x_to_libra_g09.py/g09_to_libra 
    
    f = open(filename,"r")
    A = f.readlines()
    f.close()

    # detect the lines including information from gamess output file
    info = detect_g09.detect(A,flag)

    # extract information from gamess output file
    label, Q, R = g09_extract_coordinates(A[info["coor_start"]:info["coor_end"]+1], flag)
    grad = g09_extract_gradient(A[info["grad_start"]:info["grad_end"]+1], flag)
    E_MO, C = g09_extract_mo(A[info["mo_start"]:info["mo_end"]+1], info["Ngbf"],active_space,flag)
    ao = g09_extract_ao_basis(A[info["ab_start"]:info["ab_end"]+1], label, R, flag)

    # Convert the KS excitation energies and the ground state potential energy into 
    # the total energies of excited states (1-electron basis)

    # *********Here, KS excitation energy -> total excitation energy***************
    nstates = len(states)
    E = MATRIX(nstates,nstates)
    for i in xrange(nstates):
        h_indx = states[i].from_orbit[0] - min_shift  # index of the hole orbital w.r.t. the lowest included in the active space orbital
        e_indx = states[i].to_orbit[0]   - min_shift  # --- same, only for the electron orbital
        EX_ene = info["tot_ene"] + E_MO.get(e_indx,e_indx) - E_MO.get(h_indx,h_indx) # excitation energy
        #EX_ene = info["tot_ene"] # for debugging MD without electron dynamics 
        E.set(i,i,EX_ene)

    # ******************************************************************************

    if flag == 1:
        print "ground states energy is",info["tot_ene"]
        print "Excitation Energy E is",E.show_matrix()
        print "********************************************"
        print "extract program ends"
        print "********************************************\n"
    
    return label, Q, R, grad, E, C, ao, info["Nele"]


