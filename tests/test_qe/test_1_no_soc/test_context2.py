#*********************************************************************************
#* Copyright (C) 2015-2016 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

##
# \file test_context.py This file demonstrates how to read QE wavefunctions using
# Context class of the Libra package. We also show how to use the created objects
#
#  --- No SOC case --- 
# Here, we also will consider only those orbitals that we really need
#

import os
import math

# Fisrt, we add the location of the library to test to the PYTHON path
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



print "\nTest 2: The easiest way to build some non-empty context is to read it from an XML file"

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
    gamma_only = ctx.get("Info/<xmlattr>/gamma_only","n")
    print ngw, nbnd, gamma_only

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
                coeff2.set(i+ngw, o, coeff.get(i, o).conjugate() )

        ngw = 2*ngw - 1


    if wfc_preprocess=="normalize" or wfc_preprocess=="restore":

        for i in xrange(norbs):
            mo_i = coeff2.col(i)
            nrm = (mo_i.H() * mo_i).get(0,0).real
            nrm = (1.0/math.sqrt(nrm))

            for pw in xrange(ngw):
                coeff2.set(pw,i,nrm*coeff2.get(pw,i))


    return coeff2


#################################

print "Reading energies"
e = read_qe_index("x.export/index.xml", [6,7] )
e.show_matrix()


coeff_1 = read_qe_wfc("x.export/wfc.1", "Kpoint.1", [6,7] )
coeff_2 = read_qe_wfc("x.export/wfc.2", "Kpoint.2", [6,7] )

nbnd = coeff_1.num_of_cols
print "The number of bands = ", nbnd
 
ovlp_1  = coeff_1.H() * coeff_1
ovlp_2  = coeff_2.H() * coeff_2
ovlp_12 = coeff_1.H() * coeff_2

print "<alp|alp> = "; ovlp_1.show_matrix()
print "<bet|bet> = "; ovlp_2.show_matrix()
print "<alp|bet> = "; ovlp_12.show_matrix()
print "\n"

#sys.exit(0)

for n in xrange(nbnd):
    print n, ovlp_1.get(n,n), ovlp_2.get(n,n), ovlp_12.get(n,n)  #  ovlp.get(2*n,2*n) + ovlp.get(2*n+1,2*n+1)





