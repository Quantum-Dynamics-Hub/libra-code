#*********************************************************************************
#* Copyright (C) 2017-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/


import cmath
import math
import os
import sys

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *



def compute_etot(ham, p, iM):

    ntraj = p.num_of_cols
    ndof = p.num_of_rows

    Epot, Ekin = 0.0, 0.0
    for traj in xrange(ntraj):
        Epot = Epot + ham.get_ham_adi(Py2Cpp_int([0,traj])).get(0,0).real

        for dof in xrange(ndof):
            Ekin = Ekin + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)

    Ekin = Ekin / float(ntraj) 
    Epot = Epot / float(ntraj) 
    Etot = Ekin + Epot
    
    return Ekin, Epot, Etot


def sample(x, mean_x, sigma_x, rnd):  
    nr, nc = x.num_of_rows, x.num_of_cols
    for i in range(nr):
        for j in range(nc):    
            x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() )




def compute_statistics(Q, idof, minx, maxx, dx, outfile):

    #=========== Compute grid for plotting ===========
    X = []
    sz = int((maxx - minx)/dx) + 1
    for i in range(sz):
        X.append(minx + i*dx)

    ntraj  = Q[0].num_of_cols
    pts = []
    for q in Q:
        for tr in xrange(ntraj):
            pts.append(q.get(idof, tr))

    dy1 = DATA(pts)
    dens = dy1.Calculate_Distribution(X)[0]

    f = open(outfile,"w")
    i = 0
    sz = len(X)
    for i in range(0, sz):
        f.write("%8.5f  %8.5f\n" % (X[i], dens[i] ) )
    f.close()


def compute_statistics2(Q, idof, tr, minx, maxx, dx, outfile):

    #=========== Compute grid for plotting ===========
    X = []
    sz = int((maxx - minx)/dx) + 1
    for i in range(sz):
        X.append(minx + i*dx)

    ntraj  = Q[0].num_of_cols
    pts = []
    for q in Q:
        pts.append(q.get(idof, tr))

    dy1 = DATA(pts)
    dens = dy1.Calculate_Distribution(X)[0]

    f = open(outfile,"w")
    i = 0
    sz = len(X)
    for i in range(0, sz):
        f.write("%8.5f  %8.5f\n" % (X[i], dens[i] ) )
    f.close()


