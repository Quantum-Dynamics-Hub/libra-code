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

    all_Etot = []
    Epot, Ekin = 0.0, 0.0
    for traj in xrange(ntraj):
        en = ham.get_ham_adi(Py2Cpp_int([0,traj])).get(0,0).real
        Epot = Epot + en

        ekin = 0.0
        for dof in xrange(ndof):
            ekin = ekin + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)

        Ekin = Ekin + ekin
        all_Etot.append(en+ekin)

    Ekin = Ekin / float(ntraj) 
    Epot = Epot / float(ntraj) 
    Etot = Ekin + Epot
    
    return Ekin, Epot, Etot, all_Etot


def lattice(q, Nx,Ny,Nz, dx, dy, dz, sx, sy, sz, rnd):

    ntraj = q.num_of_cols

    for traj in xrange(ntraj):
        i = 0
        for nx in xrange(Nx):
            for ny in xrange(Ny):
                for nz in xrange(Nz):

                    q.set(3*i,   traj,  nx*dx + sx * rnd.normal() )
                    q.set(3*i+1, traj,  ny*dy + sy * rnd.normal() )
                    q.set(3*i+2, traj,  nz*dz + sz * rnd.normal() )

                    i = i + 1


def lattice_p(p, Nx,Ny,Nz, px, py, pz, spx, spy, spz, rnd):

    ntraj = p.num_of_cols   
    nat = Nx * Ny * Nz
    for traj in xrange(ntraj):

        Px, Py, Pz = 0.0, 0.0, 0.0

        i = 0
        for nx in xrange(Nx):
            for ny in xrange(Ny):
                for nz in xrange(Nz):

                    p.set(3*i,   traj,  spx * rnd.normal())
                    p.set(3*i+1, traj,  spy * rnd.normal())
                    p.set(3*i+2, traj,  spz * rnd.normal())

                    Px = Px + p.get(3*i,   traj)
                    Py = Py + p.get(3*i+1, traj)
                    Pz = Pz + p.get(3*i+2, traj)
                    i = i + 1

        Px = Px/float(nat); Py = Py/float(nat); Pz = Pz/float(nat)

        i = 0
        for nx in xrange(Nx):
            for ny in xrange(Ny):
                for nz in xrange(Nz):

                    p.add(3*i,   traj,  px-Px)
                    p.add(3*i+1, traj,  py-Py)
                    p.add(3*i+2, traj,  pz-Pz)                
                    i = i + 1




def sample(x, mean_x, sigma_x, rnd):  
    nr, nc = x.num_of_rows, x.num_of_cols
    for i in range(nr):
        for j in range(nc):    
            x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() )

def chain_x(q, dx, sigma_x, sigma_y, sigma_z, rnd):  
    ndof, ntraj = q.num_of_rows, q.num_of_cols
    nat = ndof / 3

    x = 0.0
    for i in range(nat):
        for j in range(ntraj):    
            q.set(3*i,  j, x + sigma_x * rnd.normal() )
            q.set(3*i+1,j,     sigma_y * rnd.normal() )
            q.set(3*i+2,j,     sigma_z * rnd.normal() )
        x = x + dx


def chain_px(p, px, py, pz, sigma_px, sigma_py, sigma_pz, rnd):  

    ndof, ntraj = p.num_of_rows, p.num_of_cols
    nat = ndof / 3

    P = MATRIX(3, ntraj)
    for i in range(nat):
        for j in range(ntraj):    
            x = sigma_px * rnd.normal() 
            y = sigma_py * rnd.normal() 
            z = sigma_pz * rnd.normal() 
            p.set(3*i,  j, x)
            p.set(3*i+1,j, y)
            p.set(3*i+2,j, z)
            P.add(0, j, x)
            P.add(1, j, y)
            P.add(2, j, z)

    P = P / float(nat)
    # Correct momenta
    for i in range(nat):
        for j in range(ntraj):    
            p.add(3*i,  j, px - P.get(0, j))
            p.add(3*i+1,j, py - P.get(1, j))
            p.add(3*i+2,j, pz - P.get(2, j))








def print_xyz(label, q, tr, filename, t):
    ndof = q.num_of_rows
    nat = ndof / 3

    Angst = 1.889725989
    
    f = open(filename, "a")
    f.write("%5i\nStep = %5i\n" % (nat, t))
    for i in xrange(nat):
        x = q.get(3*i, tr)/Angst
        y = q.get(3*i+1, tr)/Angst
        z = q.get(3*i+2, tr)/Angst
        f.write("%s  %8.5f %8.5f %8.5f\n" % (label[i], x,y,z) )
    f.close()


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


def compute_rmsd(Q, dt, outfile):

    sz = len(Q)
    ndof = Q[0].num_of_rows
    nat = ndof/3
    ntraj = Q[0].num_of_cols


    f = open(outfile,"w")
    for i in xrange(sz):

        rmsd = 0.0
        for tr in xrange(ntraj):
            dQ = Q[i].col(tr) - Q[0].col(tr) 
            rmsd = rmsd + (dQ.T() * dQ).get(0)

        rmsd = rmsd / float(ntraj * nat)

        f.write("%8.5f  %8.5f\n" % (i*dt, rmsd ) )
    f.close()

    
def compute_cv(e, start, dt, outfile):
    """
    e - [ [E(t=0)], [E(t=1)], ...],
        each E(t) is [ [E(traj=0)], [E(traj=1)], ... ]

    start - how many initial datapoints to discard
    dt - integration time step
    outfile - where to print out the data
    """

    sz = len(e)
    ntraj = len(e[0])

    ave = MATRIX(1, ntraj)
    for i in xrange(sz):
        for traj in xrange(ntraj):
            ave.add(0, traj, e[i][traj])
    ave = ave/float(sz)
   

    # Compute fluctuation
    de2 = MATRIX(sz, 1)
    for i in xrange(sz):        
        tot = 0.0
        for traj in xrange(ntraj):
            de  = e[i][traj] - ave.get(0, traj)
            tot = tot + de * de
        tot = tot/float(ntraj)
        de2.set(i, 0, tot)
     
    f = open(outfile,"w")
    f.close()

    f = open(outfile,"a")
    cum = 0.0
    for i in xrange(sz):
        if i>start:
            n = i - start
            cum = (n*cum + de2.get(i) )/float(n+1.0)
        else:
            cum = de2.get(i)

        f.write("%8.5f %8.5f %8.5f \n" % (i*dt, de2.get(i), cum ))
    f.close()

    return cum

