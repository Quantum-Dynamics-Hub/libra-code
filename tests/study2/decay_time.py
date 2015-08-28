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

import math

def Regression(X,Y,sz):
# Finds y = a + b*x
    x,y,x2,y2,xy = 0.0, 0.0, 0.0, 0.0, 0.0

    i = 0
    while i<sz:
        x = x + X[i]
        y = y + Y[i]
        x2 = x2 + X[i]*X[i]
        y2 = y2 + Y[i]*Y[i]
        xy = xy + X[i]*Y[i]
        i = i + 1
    N = float(sz)
    b = (N*xy - x*y)/(N*x2 - x*x)
    a = (y - b*x)/N
#    b = y/x
#    a = 0.0

    return [a,b]

 
for indx in range(0,17):

    filename = "relax"+str(indx)+".txt"

    f = open(filename,"r")
    A = f.readlines()
    f.close()

    T, P = [], []
    sz = len(A)
    act_sz = 0
    i = 0
    while i<sz:
        tmp = A[i].split()

        t = float(tmp[0])
        val = float(tmp[1])

        if val>0.01: #and t<2000000:
            T.append(t)
            P.append(math.log(val) ) 
            act_sz = act_sz + 1
        i = i + 1

    a,b = Regression(T,P,act_sz)

    tau = -1.0/b
    A = math.exp(a)
    print "tau = ",tau, "k = ", -b, "A = ", A

    f = open("model"+str(indx),"w")

    for t in T:
        p = A * math.exp(-t/tau)
        line = str(t)+"   "+str(p)+"\n"
        f.write(line)
    f.close()




    