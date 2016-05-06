#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
print os.getcwd()
cwd = "/home/Alexey_2/Programming/Project_libra/_build" #os.getcwd()
print "Current working directory", cwd
sys.path.insert(1,cwd+"/src/mmath")
sys.path.insert(1,cwd+"/src/context")


print "\nTest 1: Importing the library and its content"
from cygmmath import *
from cygcontext import *


print "\nTest 2: The easiest way to build some non-empty context is to read it from an XML file"
ctx = Context("x.export/wfc.1")
ctx.set_path_separator("/")
print "path=", ctx.get_path()
ctx.save_xml("x.export/wfc1.xml")
ctx.show_children("Kpoint.1")

#sys.exit(0)

ngw = int(float(ctx.get("Info/<xmlattr>/ngw","n")))
nbnd = int(float(ctx.get("Info/<xmlattr>/nbnd","n")))

coeff = CMATRIX(ngw,nbnd)

#print "Verbatim string \n", ctx.get("Wfc.1", "n")

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


#print "Matrix .."
#coeff.show_matrix()

ovlp = CMATRIX(nbnd, nbnd)
ovlp = coeff.H() * coeff

#ovlp.show_matrix()

for n in xrange(nbnd/2):
    print n, ovlp.get(2*n,2*n) + ovlp.get(2*n+1,2*n+1)

for n in xrange(nbnd):
    print n, ovlp.get(n,n)





