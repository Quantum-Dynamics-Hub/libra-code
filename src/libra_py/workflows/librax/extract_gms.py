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

## \file extract_gms.py
# This module implements the functions that extract
# atomic forces , molecular energies, molecular orbitals, and atomic basis information
# written in gamess output file.

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

import detect_gms
import ao_basis

def gms_extract_ao_basis(inp_str, label, R, flag):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] inp_str A list containing info. for atomic orbital basis.
    # \param[in] label   A list of atomic labels (e.g. C,H,O)
    # \param[in] R       A list of atomic coordinates (R[i] =[R[i].x, R[i].y, R[i].z]
    # \param[in] flag    debugging info: option 1-> print, otherwise -> don't 
    # ao - returned lists of atomic orbital basis sets 
    #
    # Used in: extract_gms.py/gms_extract

    # atomic species
    l_atom_spec = []
    atom_spec = []
    sz = len(inp_str)

    # atom labels
    for i in xrange(sz):
        spline = inp_str[i].split()
        if len(spline) == 1:
            l_atom_spec.append(i)
            # **** Here, atomic names are changed *****
            # atom labels
            atom = spline[0]
            if len(atom) > 1:
                atom = atom[0] + atom[1].lower()
            atom_spec.append(atom)
            # *****************************************


    # atomic basis sets
    basis_type = []
    basis_expo = []
    basis_coef = []

    for i in range(0,len(atom_spec)):
        stmp = atom_spec[i]
        type_tmp = []
        expo_tmp = []
        coef_tmp = []
        if i < len(atom_spec) - 1:
            tmp_start = l_atom_spec[i] + 2
            tmp_end = l_atom_spec[i+1] - 2
        else:
            tmp_start = l_atom_spec[i] + 2
            tmp_end = sz - 1 #ab_end

        for j in range(tmp_start,tmp_end+1):
            spline = inp_str[j].split()
            coef_tmp1 = []
            if len(spline) > 1:
                type_tmp.append(spline[1])
                expo_tmp.append(float(spline[3]))
                if spline[1] == "L":
                    coef_tmp1.append(float(spline[4]))
                    coef_tmp1.append(float(spline[5]))
                else:
                    coef_tmp1.append(float(spline[4]))
                coef_tmp.append(coef_tmp1)
        basis_type.append(type_tmp)
        basis_expo.append(expo_tmp)
        basis_coef.append(coef_tmp)


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
            if b_tmp == "S":
                expo_stmp.append(basis_expo[i][j])
                coef_stmp.append(basis_coef[i][j][0])
            elif b_tmp == "P":
                expo_ptmp.append(basis_expo[i][j])
                coef_ptmp.append(basis_coef[i][j][0])
            elif b_tmp == "D":
                expo_dtmp.append(basis_expo[i][j])
                coef_dtmp.append(basis_coef[i][j][0])
            elif b_tmp == "L":
                expo_stmp.append(basis_expo[i][j])
                expo_ptmp.append(basis_expo[i][j])
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


def gms_extract_mo(inp_str,Ngbf,active_space,flag):
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
    # Used in: extract_gms.py/gms_extract

    stat_span = Ngbf + 4 # period for beginning of coefficient lines

    sz = len(inp_str) 

    # create objects of MATRIX type, containing eigenvalues and eigenvectors 
    E_full = MATRIX(Ngbf,Ngbf)
    C_full = MATRIX(Ngbf,Ngbf)

    for i in xrange(sz):

        if i % stat_span == 0 :
            ind_of_eig = inp_str[i].split() # split lines for indexes of eigenvalues

            # eigenvalues
            eig_val = inp_str[i+1].split() # split lines for eigenvalues
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
                    if len(eig_vec[1]) == 4:        # continuous word like "CL44"
                        kvec = k + 3
                    else:                           # discontinuous word like "H 20" 
                        kvec = k + 4
                    C_full.set(ja,ke,float(eig_vec[kvec]))

    #Here, we need to add reduction of E, C matrices and the corresponding trunctaion of the ao list. 
    #This will allow us make further computations faster and in consistent with those adopted in QE 

    # ***********Here, reduce E_full and C_full ***************
    sz = len(active_space)
    if sz==0:
        print "active space is not defined correctly, exit....."
        sys.exit(0)

    #E = MATRIX(sz,sz)
    #for i in xrange(sz):
    #    imo = active_space[i]-1
    #    E.set(i,i,E_full.get(imo,imo))

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
        #print "E Matrix is"; E.show_matrix()
        print "C Matrix is"; C.show_matrix()

    return E_full, C
    
def gms_extract_coordinates(inp_str,flag):
    ##
    # Extracts atomic labels, nuclear charges, and coordinates of all atoms
    # from the the list of input lines
    # each input line is assumed to have the format:
    # label  Q  ....  x y z
    #
    # \param[in] inp_str  Strings containing the info for all atomic coordinates
    # \param[in] flag     debugging info: option 1-> print, otherwise -> don't
    # label - returned list of atomic labels (strings)
    # Q - returned list of nuclear charges (floats)
    # R - returned list of nuclear coordinates (VECTOR objects)
    #
    # Used in: extract_gms.py/gms_extract


    label, Q, R = [], [], []
    for a in inp_str: 
        spline = a.split() 

        # atom labels
        atom = spline[0]
        if len(atom) > 1:
            atom = atom[0] + atom[1].lower()
        label.append(atom)

        # atomic charges
        Q.append(float(spline[1]))

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


def gms_extract_gradient(inp_str,flag):
    ##
    # Extracts atomic gradients on all atoms from the the list of input lines
    # each input line is assumed to have the format:
    #  ....  gx  gy  gz
    #
    # \param[in] inp_str  Strings containing the gradient for all atoms
    # \param[in] flag     debugging info: option 1-> print, otherwise -> don't
    # grad - returned list of VECTOR objects
    #
    # Used in: extract_gms.py/gms_extract

    grad = []
    for a in inp_str:
        spline = a.split()  

        # Last 3 elements of the line
        x = float(spline[len(spline)-3])
        y = float(spline[len(spline)-2])
        z = float(spline[len(spline)-1])
        g = VECTOR(x,y,z)

        grad.append(g)

    if flag == 1:
        print "atomic gradient="
        for g in grad:
            print grad.index(g), g, g.x, g.y, g.z

    return grad


def gms_extract(filename,states,min_shift,active_space,flag):
    ##
    # This function extracts all the necessary information (energies, gradients, coordinates, MOs, 
    # AOs, etc. ) from the GAMESS output.
    #
    # \param[in] filename    The name of the GAMESS output file from which we will be getting info
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
    # Used in: main.py/main or x_to_libra_gms.py/gamess_to_libra

    f = open(filename,"r")
    A = f.readlines()
    f.close()

    # detect the lines including information from gamess output file
    info = detect_gms.detect(A,flag)

    # extract information from gamess output file
    label, Q, R = gms_extract_coordinates(A[info["coor_start"]:info["coor_end"]+1], flag)
    grad = gms_extract_gradient(A[info["grad_start"]:info["grad_end"]+1], flag)
    E_MO, C = gms_extract_mo(A[info["mo_start"]:info["mo_end"]+1], info["Ngbf"],active_space,flag)
    ao = gms_extract_ao_basis(A[info["ab_start"]:info["ab_end"]+1], label, R, flag)

    # Convert the KS excitation energies and the ground state potential energy into 
    # the total energies of excited states (1-electron basis)

    # *********Here, KS excitation energy -> total excitation energy***************
    homo = active_space[0] - min_shift # HOMO : number ordering starts from 1, not 0.
    nstates = len(states)
    E = MATRIX(nstates,nstates)
    for i in xrange(nstates):
        #h_indx = states[i].from_orbit[0] - min_shift  # index of the hole orbital w.r.t. the lowest included in the active space orbital
        #e_indx = states[i].to_orbit[0]   - min_shift  # --- same, only for the electron orbital
        #print "%i th excitation" %(i)
        h_indx = states[i].from_orbit[0] + homo - 1
        e_indx = states[i].to_orbit[0] + homo - 1
        EX_ene = info["tot_ene"] + E_MO.get(e_indx,e_indx) - E_MO.get(h_indx,h_indx) # excitation energy
        #EX_ene = info["tot_ene"] # for debugging MD without electron dynamics 
        E.set(i,i,EX_ene)
        if 0==1: # debug
            print "e_indx is %i" %(e_indx)
            print "MO energy is %f" %( E_MO.get(e_indx,e_indx))
            print "h_indx is %i" %(h_indx)
            print "MO energy is %f" %( E_MO.get(h_indx,h_indx))
    # ******************************************************************************

    if flag == 1:
        print "ground states energy is",info["tot_ene"]
        print "Excitation Energy E is",E.show_matrix()
        print "********************************************"
        print "extract program ends"
        print "********************************************\n"
    
    return label, Q, R, grad, E, C, ao, info["Nele"]


