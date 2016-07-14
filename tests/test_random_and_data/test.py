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

import os
import math

# Fisrt, we add the location of the library to test to the PYTHON path
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


print "\nTest 1: Constructor is called each time when the number is generated"
for i in xrange(10):
    r = Random()
    print i, r.uniform(0.0, 1.0)

print "\nTest 1a: Constructor is created only once"
r = Random()
for i in xrange(10):
    print i, r.uniform(0.0, 1.0)



print "\nTest 2: Constructor & uniform"
r = Random()

y1 = []
for i in range(0,100000):
    y1.append( r.uniform(0.0, 1.0) )
dy1 = DATA(y1)

x = []
for i in range(0,100):
    x.append(i*0.01)

dens, cum = dy1.Calculate_Distribution(x)

f = open("uniform.txt","w")
i = 0
sz = len(x)
for i in range(0, sz):
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  \n" % (x[i], dens[i], cum[i], r.p_uniform(0.0, 1.0)) )
f.close()



print "\nTest 3: exponential"
y1 = []
for i in range(0,100000):
    y1.append( r.exponential(1.0) )
dy1 = DATA(y1)

x = []
for i in range(0,100):
    x.append(i*0.1)

dens, cum = dy1.Calculate_Distribution(x)

f = open("exponential.txt","w")
i = 0
sz = len(x)
for i in range(0, sz):
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  \n" % (x[i], dens[i], cum[i], r.p_exponential(x[i], 1.0)) )
f.close()




print "\nTest 4: normal"
y1 = []
for i in range(0,100000):
    y1.append( r.normal() )
dy1 = DATA(y1)

x = []
for i in range(-100,100):
    x.append(i*0.1)

dens, cum = dy1.Calculate_Distribution(x)

f = open("normal.txt","w")
i = 0
sz = len(x)
for i in range(0, sz):
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  \n" % (x[i], dens[i], cum[i], r.p_normal(x[i])) )
f.close()



print "BEWARE!!!  Gamma and beta functions need fixing! The distributions are probably right, but the PDF are not as they should be\n";
print "\nTest 5: gamma"
y1 = []
for i in range(0,100000):
    y1.append( r.gamma(6) )
dy1 = DATA(y1)

x = []
for i in range(0,100):
    x.append(i*0.1)

dens, cum = dy1.Calculate_Distribution(x)

f = open("gamma.txt","w")
i = 0
sz = len(x)
for i in range(0, sz):
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  \n" % (x[i], dens[i], cum[i], r.p_gamma(1.2, x[i])) )
f.close()



print "\nTest 6: poisson (poiss1)"
y1 = []
for i in range(0,1000000):
    y1.append( r.poiss1(5.0) )
dy1 = DATA(y1)

x = []
for i in range(0,25):
    x.append(i)

dens, cum = dy1.Calculate_Distribution(x)

f = open("poiss1.txt","w")
i = 0
sz = len(x)
for i in range(0, sz):
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  \n" % (x[i], dens[i], cum[i], r.p_poiss(x[i], 5.0)) )
f.close()



print "\nTest 6a: poisson (poiss2)"
y1 = []
for i in range(0,1000000):
    y1.append( r.poiss2(1.0) )
dy1 = DATA(y1)

x = []
for i in range(0,10):
    x.append(i)

dens, cum = dy1.Calculate_Distribution(x)

f = open("poiss2.txt","w")
i = 0
sz = len(x)
for i in range(0, sz):
    f.write("%8.5f  %8.5f  %8.5f  %8.5f  \n" % (x[i], dens[i], cum[i], r.p_poiss(x[i], 1.0)) )
f.close()





