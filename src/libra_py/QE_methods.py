#*********************************************************************************
#* Copyright (C) 2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file QE_methods.py
# This module implements functions for dealing with the outputs from QE (Quantum Espresso) package

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def read_qe_index(filename, orb_list, verbose=0):
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the coefficients of the plane waves that constitute the
#  wavefunction
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#  \param[in] upper_tag This is the name of the upper-level tag
#  \param[in] orb_list The list containing the indices of the orbitals which we want to consider
#   (indexing is starting with 1, not 0!)
#  
# Returns: a diagonal matrix (complex-valued) of orbital energies - only for those that are of interest
# in Ha = a.u. of energy


    Ry2Ha = 0.5 # conversion factor
   
    ctx = Context(filename)  #("x.export/index.xml")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()
    ctx.show_children("Root/Eigenvalues")  #("Kpoint.1")

    info = {}
    info["nspin"] = int(float(ctx.get("Eigenvalues/<xmlattr>/nspin","-1.0")))  # 1 - non-polarized, 2 - polarized, 4 - non-collinear
    info["nk"] = int(float(ctx.get("Eigenvalues/<xmlattr>/nk","-1.0")))  # the number of k-points
    info["nbnd"] = int(float(ctx.get("Eigenvalues/<xmlattr>/nbnd","-1.0")))  # the number of orbitals in each k-point
    info["efermi"] = float(ctx.get("Eigenvalues/<xmlattr>/efermi","-1.0"))  # Fermi energy

    # Usially in atomic units!
    info["alat"] = float(ctx.get("Cell/Data/<xmlattr>/alat","-1.0"))    # lattice constant (a)
    info["omega"] = float(ctx.get("Cell/Data/<xmlattr>/omega","-1.0"))  # unit cell volume
    info["tpiba"] = float(ctx.get("Cell/Data/<xmlattr>/tpiba","-1.0"))  # 2 pi / a
    info["tpiba2"] = float(ctx.get("Cell/Data/<xmlattr>/tpiba2","-1.0"))  # ( 2 pi / a )**2

    # Direct lattice vectors
    tmp = ctx.get("Cell/a1/<xmlattr>/xyz","-1.0").split()
    info["a1"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/a2/<xmlattr>/xyz","-1.0").split()
    info["a2"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/a3/<xmlattr>/xyz","-1.0").split()
    info["a3"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))


    # Reciprocal lattice vectors
    tmp = ctx.get("Cell/b1/<xmlattr>/xyz","-1.0").split()
    info["b1"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/b2/<xmlattr>/xyz","-1.0").split()
    info["b2"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))

    tmp = ctx.get("Cell/b3/<xmlattr>/xyz","-1.0").split()
    info["b3"] = VECTOR(float(tmp[0]), float(tmp[1]), float(tmp[2]))


    # K-points
    # Weights of the k-points
    tmp = ctx.get("Kmesh/weights","-1.0").split()
    weights = []
    for x in tmp:
        weights.append( float(x) )
    info["weights"] = weights

    # K-point vectors
    K = []
    tmp = ctx.get("Kmesh/k","-1.0").split()
    nk = len(tmp)/3 # the number of k-points
    for ik in xrange(nk):
        k = VECTOR(float(tmp[3*ik+0]), float(tmp[3*ik+1]), float(tmp[3*ik+2]))
        K.append(k)
    info["k"] = K       

    
    # All the bands
    all_e = []  # all bands

    for ik in range(1,info["nk"]+1):

        e_str = ctx.get("Eigenvalues/e.%i" % ik,"n").split() 
        nbnd = len(e_str)    

        norbs = len(orb_list)
        for o in orb_list:
            if o > nbnd:
                print "Orbital ", o, " is outside the range of allowed orbital indices. The maximal value is ", nbnd
                sys.exit(0)
            elif o < 1:
                print "Orbital ", o, " is outside the range of allowed orbital indices. The minimal value is ", 1
                sys.exit(0)


        e = CMATRIX(norbs,norbs)
        
        for band in orb_list:    
            mo_indx = orb_list.index(band) 
        
            ei = float(e_str[band-1])  # band starts from 1, not 0!
        
            e.set(mo_indx,mo_indx, ei * Ry2Ha, 0.0)

        all_e.append(e)

    if verbose==1:
        print " nspin = ", info["nspin"], " nk = ", info["nk"], " nbnd = ", info["nbnd"], " efermi = ", info["efermi"], \
        " alat = ", info["alat"], " omega = ", info["omega"], " tpiba = ", info["tpiba"], " tpiba2 = ", info["tpiba"]

        print " Direct lattice vectors: "
        print " a1 = ", info["a1"].x, info["a1"].y, info["a1"].z
        print " a2 = ", info["a2"].x, info["a2"].y, info["a2"].z
        print " a3 = ", info["a3"].x, info["a3"].y, info["a3"].z

        print " Reciprocal lattice vectors: "
        print " b1 = ", info["b1"].x, info["b1"].y, info["b1"].z
        print " b2 = ", info["b2"].x, info["b2"].y, info["b2"].z
        print " b3 = ", info["b3"].x, info["b3"].y, info["b3"].z

        print " K points: "
        for ik in xrange(info["nk"]):
            print ik, " weight = ", info["weights"][ik], " k = ", info["k"][ik].x, info["k"][ik].y, info["k"][ik].z

        print " Energies of the active orbitals for all k-points: "
        for ik in xrange(info["nk"]):
            print "ik = ", ik
            all_e[ik].show_matrix()
            print ""


    return info, all_e



def read_qe_wfc_info(filename, verbose=0):
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the some descriptors
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#  
   
    ctx = Context(filename)  #("x.export/wfc.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()

    res = {}
    res["ngw"] = int(float(ctx.get("Info/<xmlattr>/ngw","-1.0")))     # the number of plane waves needed to represent the orbital
    res["igwx"] = int(float(ctx.get("Info/<xmlattr>/igwx","-1.0")))   # the number of the G points = plane waves needed to
                                                                      # represent the orbital for given k-point. Use this number 
                                                                      # when working with multiple k-points
    res["nbnd"] = int(float(ctx.get("Info/<xmlattr>/nbnd","-1.0")))   # the number of bands (orbitals)
    res["nspin"] = int(float(ctx.get("Info/<xmlattr>/nspin","-1.0"))) # 1 - unpolarized, 2 - polarized, 4 - non-collinear 
    res["gamma_only"] = ctx.get("Info/<xmlattr>/gamma_only","F")      # T - use the Gamma-point storae trick, T - do not use it
    res["ik"] = int(float(ctx.get("Info/<xmlattr>/ik","-1.0")))       # index of the k point wfc
    res["nk"] = int(float(ctx.get("Info/<xmlattr>/nk","-1.0")))       # the number of k-points in the wfc


    if verbose==1:
        print " ngw = ", res["ngw"], " igwx = ", res["igwx"],\
              " nbnd = ", res["nbnd"], " nspin = ", res["nspin"],\
              " gamma_only = ", res["gamma_only"],\
              " ik = ", res["ik"], " nk = ", res["nk"]

    return res




def read_qe_wfc_grid(filename, verbose=0):
##
#  This functions reads an ASCII/XML format file containing grid of G-points for given k-point
#  and returns a list of VECTOR objects
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#     

    ctx = Context(filename)  #("x.export/grid.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()


    # G-point vectors
    G = []
    tmp = ctx.get("grid","-1.0").split()
    ng = len(tmp)/3 # the number of G-points

    for ig in xrange(ng):
        g = VECTOR(float(tmp[3*ig+0]), float(tmp[3*ig+1]), float(tmp[3*ig+2]))
        G.append(g)

    if verbose==1:
        print "The number of G point on the grid for this k-point = ", ng
        print " G points: "
        for ig in xrange(ng):
            print ig, " g = ", G[ig].x, G[ig].y, G[ig].z
    

    return G



def read_qe_wfc(filename, orb_list, verbose=0):
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the coefficients of the plane waves that constitute the
#  wavefunction
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#  \param[in] orb_list The list containing the indices of the orbitals which we want to consider
#   (indexing is starting with 1, not 0!)
#     

    ctx = Context(filename)  #("x.export/wfc.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()

    res = read_qe_wfc_info(filename, verbose)
    ngw, nbnd, nspin, gamma_only = res["igwx"], res["nbnd"], res["nspin"], res["gamma_only"]

    if nspin==4:
        if orb_list[0] % 2 ==0:
            print "In SOC, the very first orbital index must be odd!\nExiting now..."
            sys.exit(0)
        if len(orb_list) % 2 == 1:
            print "In SOC, an even number of orbitals must be utilized!\nExiting now..."
            sys.exit(0)

    wfc_preprocess = "normalize"
    if gamma_only=="T":
        wfc_preprocess = "restore"


    norbs = len(orb_list)
    for o in orb_list:
        if o > nbnd:
            print "Orbital ", o, " is outside the range of allowed orbital indices. The maximal value is ", nbnd
            sys.exit(0)
        elif o < 1:
            print "Orbital ", o, " is outside the range of allowed orbital indices. The minimal value is ", 1
            sys.exit(0)


    coeff = CMATRIX(ngw,norbs)

    for band in orb_list:
        c = []
        all_coeff = ctx.get("Wfc."+str(band), "n").split(',')
        sz = len(all_coeff)

        for i in xrange(sz):
            a = all_coeff[i].split()
            for j in xrange(len(a)):
                c.append(a[j])
        sz = len(c)
        n = sz/2  # this should be equal to ngw

        mo_indx = orb_list.index(band) # - 1
        for i in xrange(n):
            coeff.set(i, mo_indx, float(c[2*i]), float(c[2*i+1]))



    #======== Normalize or restore (gamma-point trick) wavefunction ===============
    coeff2 = coeff #CMATRIX(ngw,norbs)

    if wfc_preprocess=="restore":

        coeff2 = CMATRIX(2*ngw-1,norbs)
        for o in xrange(norbs):
            coeff2.set(0, o, coeff.get(0, o))

            for i in xrange(1,ngw):
                coeff2.set(i, o, coeff.get(i, o))
                coeff2.set(i+ngw-1, o, coeff.get(i, o).conjugate() )

        ngw = 2*ngw - 1


    if wfc_preprocess=="normalize" or wfc_preprocess=="restore":

        if nspin==1 or nspin==2:

            for i in xrange(norbs):
                mo_i = coeff2.col(i)            
                nrm = (mo_i.H() * mo_i).get(0,0).real
                nrm = (1.0/math.sqrt(nrm))
        
                for pw in xrange(ngw):
                    coeff2.set(pw,i,nrm*coeff2.get(pw,i))

        elif nspin==4:  # spinor case

            for i in xrange(norbs/2):  # this will always be even
                mo_i_a = coeff2.col(2*i)
                mo_i_b = coeff2.col(2*i+1)

                nrm = ( (mo_i_a.H() * mo_i_a).get(0,0).real + (mo_i_b.H() * mo_i_b).get(0,0).real )
                nrm = (1.0/math.sqrt(nrm))
        
                for pw in xrange(ngw):
                    coeff2.set(pw,2*i,nrm*coeff2.get(pw,2*i))
                    coeff2.set(pw,2*i+1,nrm*coeff2.get(pw,2*i+1))
        

    return coeff2


