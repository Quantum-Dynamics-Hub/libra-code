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
#  --- SOC case --- 

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

def read_qe_wfc(filename, upper_tag):
##
#  This functions reads an ASCII/XML format file containing wavefunction
#  and returns the coefficients of the plane waves that constitute the
#  wavefunction
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
    print ngw, nbnd

    coeff = CMATRIX(ngw,nbnd)

    for band in range(1,nbnd+1):
        c = []
        all_coeff = ctx.get("Wfc."+str(band), "n").split(',')
        sz = len(all_coeff)

        for i in xrange(sz):
            a = all_coeff[i].split()
            for j in xrange(len(a)):
                c.append(a[j])
        sz = len(c)
        n = sz/2  # this should be equal to ngw

        for i in xrange(n):
            coeff.set(i, band-1, float(c[2*i]), float(c[2*i+1]))

    return coeff

#################################


coeff_1 = read_qe_wfc("x.export/wfc.1", "Kpoint.1")

nbnd = coeff_1.num_of_cols
print "The number of bands = ", nbnd
 
ovlp_1  = coeff_1.H() * coeff_1

print "<mix|mix> = "; ovlp_1.show_matrix()
print "\n"

for n in xrange(nbnd):
    print n, ovlp_1.get(n,n)
print "\n"

for n in xrange(nbnd/2):
    print n, ovlp_1.get(2*n,2*n) + ovlp_1.get(2*n+1,2*n+1)
print "\n"

ovlp = CMATRIX(nbnd/2, nbnd/2)
for n in xrange(nbnd/2):
    for k in xrange(nbnd/2):
        ovlp.set(n,k, ovlp_1.get(2*n,2*k) + ovlp_1.get(2*n+1,2*k+1) )

print "spinor overlap matrix = "; ovlp.show_matrix()


