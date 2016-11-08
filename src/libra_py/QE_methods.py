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


def read_qe_index(filename, orb_list):
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


    e_str = ctx.get("Eigenvalues/e.1","n").split() 
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


    return e



def read_qe_wfc_info(filename, upper_tag):
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the some descriptors
#
#  \param[in] filename This is the name of the file we will be reading to construct a wavefunction
#  \param[in] upper_tag This is the name of the upper-level tag
#  
   
    ctx = Context(filename)  #("x.export/wfc.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()
    ctx.show_children(upper_tag)  #("Kpoint.1")

    ngw = int(float(ctx.get("Info/<xmlattr>/ngw","n")))
    nbnd = int(float(ctx.get("Info/<xmlattr>/nbnd","n")))
    nspin = int(float(ctx.get("Info/<xmlattr>/nspin","n")))
    gamma_only = ctx.get("Info/<xmlattr>/gamma_only","n")

    print "ngw = ", ngw, " nbnd = ", nbnd, " nspin = ", nspin, "gamma_only = ", gamma_only


    return ngw, nbnd, nspin, gamma_only




def read_qe_wfc(filename, upper_tag, orb_list):
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
   
    ctx = Context(filename)  #("x.export/wfc.1")
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()
    ctx.show_children(upper_tag)  #("Kpoint.1")

    ngw = int(float(ctx.get("Info/<xmlattr>/ngw","n")))
    nbnd = int(float(ctx.get("Info/<xmlattr>/nbnd","n")))
    nspin = int(float(ctx.get("Info/<xmlattr>/nspin","n")))
    gamma_only = ctx.get("Info/<xmlattr>/gamma_only","n")
    print "ngw = ", ngw, " nbnd = ", nbnd, " nspin = ", nspin, "gamma_only = ", gamma_only

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


