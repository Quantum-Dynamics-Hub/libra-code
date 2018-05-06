#*********************************************************************************
#* Copyright (C) 2017-2018  Brendan A. Smith, Alexey V. Akimov
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


def check_potential(q, params, minx_, maxx_, dx, miny_, maxy_, dy, filename):
    """
    This specific example generates a 2D potential energy surface, and checks
    if Martens' Model 2 2D potential is impimented correctly.

    q - matrix of trajector positions, of sinze ndof x ntraj

    Martens' model I  = Va*sech(2.0*x)*sech(2.0*x) + 0.5*Vb*y*y
    Martens' model II = Va*sech(2.0*x)*sech(2.0*x) + 0.5*Vb*(y + Vc*(x*x-1.0))**2)
    """
    
    # Prepare the grids, for each dimension    
    npts_x = int((maxx_ - minx_)/dx) + 1 
    npts_y = int((maxy_ - miny_)/dy) + 1 

    range_x = list(i*0.1 for i in xrange(-npts_x/2,npts_x/2))
    range_y = list(i*0.1 for i in xrange(-npts_y/2,npts_y/2))

    # Define parameters
    model = params["model"]
    Va = 0.00625
    Vb = 0.0106
    Vc = 0.4

    pes = MATRIX(npts_x,npts_y)

    for i in xrange(npts_x) :
        for j in xrange(npts_y):

            x = range_x[i]
            y = range_y[j]  

            x2 = x*x
            y2 = y*y

            if model == 1:
                pot = 0.25*(x2*x2 + y2*y2) - 0.5*(x2 + y2)
            if model == 2: 
                pot = Va*sech(2.0*x)*sech(2.0*x) + 0.5*Vb*(y+Vc*(x2 - 1.0))**2
            if model == 3:
                pot = Va*sech(2.0*x)*sech(2.0*x) + 0.5*Vb*y2

            pes.add(i,j,pot)

    datautils.show_matrix_splot(pes,filename)    

def sech(x):

    res = 1.0/math.cosh(x)

    return res

def compute_etot(ham, p, iM):

    ntraj = p.num_of_cols
    ndof = p.num_of_rows

    Epot, Ekin = 0.0, 0.0
    for traj in xrange(ntraj):
        Epot = Epot + ham.get_ham_adi(Py2Cpp_int([0,traj])).get(0,0).real
        for dof in xrange(ndof):          
            Ekin = Ekin + 0.5 * iM.get(dof, 0) * (p.get(dof, traj) ** 2)

    for i in xrange(ndof):
        Epot += ham.get_ham_adi().get(i,0).real

    Epot /=  (float(ntraj))
    Ekin /=  (float(ntraj))

    Etot = Ekin + Epot
    
    return Ekin, Epot, Etot

def sample(x, mean_x, sigma_x, rnd, sample_opt):  
    nr, nc = x.num_of_rows, x.num_of_cols

    if sample_opt == 0:
        for i in range(nr):
            for j in range(nc):    
                x.set(i,j, mean_x.get(i,0) + sigma_x.get(i,0) * rnd.normal() )

    elif sample_opt == 1:
        params = {"k":0.032, "m":2000.0, "states":[0], "coeffs":[1.0]}
        sampling = metropolis_gau(rnd, HO_sup, x, params, 1000, 500, 0.25)

        for i in range(nr):
            for j in range(nc):
                    x.set(i,j, sampling[j].get(0,i) + mean_x.get(i,0))

    elif sample_opt == 2:
        params = {"k":0.032, "m":2000.0, "states":[0], "coeffs":[1.0]}
        sampling = metropolis_gau(rnd, HO_sup, x, params, 1000, 500, 2.0)

        for i in range(nr):
            for j in range(nc):
                    x.set(i,j, sampling[j].get(0,i) + mean_x.get(i,0))




def bin(sample, min_, max_, dx, filename):
    """
    sample = dimentions are: (nnucl x ntraj)
    """
    ndof, ntraj = sample.num_of_rows, sample.num_of_cols

    # Prepare the grids, for each dimension
    max_pts = int((max_ - min_)/dx) + 1
    x_points, y_points = MATRIX(ndof,max_pts), MATRIX(ndof,max_pts)
    for i in xrange(ndof):
        for j in xrange(max_pts):
            x_points.set(i,j,min_ + j * dx)
            y_points.set(i,j,0.0)

    # Compute the frequencies 
    for i in xrange(ndof):
        for j in xrange(ntraj):
            x = sample.get(i,j) 
            indx = int((x - min_)/dx)
            y_points.set( i,indx,(y_points.get(i,indx) + 1.0/float(ntraj)) )
 
    for i in xrange(ndof):
        f = open(filename+"_distrib_dof_"+str(i)+".txt", "w")     
        for j in xrange(max_pts):
            f.write("%8.5f  %8.5f \n" % (x_points.get(i,j), y_points.get(i,j)) )
    f.close()


def bin2(sample, minx_, maxx_, dx, miny_, maxy_, dy, filename):
    """
    sample = dimentions are: (nnucl x ntraj)
    """

    ndof, ntraj = sample.num_of_rows, sample.num_of_cols

    # Prepare the grids, for each dimension
    npts_x = int((maxx_ - minx_)/dx) + 1
    npts_y = int((maxy_ - miny_)/dy) + 1

    distr = MATRIX(npts_x, npts_y) 

    for traj in xrange(ntraj):      
   
        indx_x = int( (sample.get(0,traj) - minx_)/dx )
        indx_y = int( (sample.get(1,traj) - miny_)/dy )

        distr.add(indx_x,indx_y, 1.0)

    distr *= (1.0/float(ntraj))

    # Printing is taken care of via Libra.py module
    datautils.show_matrix_splot(distr,filename)


def extract_q_p_info(q,p):

    ntraj = q.num_of_cols

    q0, q1 = [], []
    p0, p1 = [], []
    for traj in xrange(ntraj):
            q0.append( q.get(0,traj) )
            q1.append( q.get(1,traj) )
            p0.append( p.get(0,traj) )
            p1.append( p.get(1,traj) )

    os.system("mkdir _q_p_info")

    out1 = open("_q_p_info/_q0_info.txt", "w"); out1.close()
    out2 = open("_q_p_info/_q1_info.txt", "w"); out2.close()
    out3 = open("_q_p_info/_p0_info.txt", "w"); out3.close()
    out4 = open("_q_p_info/_p1_info.txt", "w"); out4.close()

    for traj in xrange(ntraj):

        out1 = open("_q_p_info/_q0_info.txt", "a")
        out2 = open("_q_p_info/_q1_info.txt", "a")
        out3 = open("_q_p_info/_p0_info.txt", "a")
        out4 = open("_q_p_info/_p1_info.txt", "a")

        out1.write( "%s\n" % ( q0[traj] ) )
        out2.write( "%s\n" % ( q1[traj] ) )
        out3.write( "%s\n" % ( p0[traj] ) )
        out4.write( "%s\n" % ( p1[traj] ) )

        out1.close()
        out2.close()
        out3.close()
        out4.close()

def get_q_p_info(params):
    """
    This function extracts the initial position and momenta coordinates
    from a previous simulation, in which the user opted to store the initial
    position and momenta coordinates

    Returns the q and p coodiates in MATRIX form
    """

    aa = open("_q_p_info/_q0_info.txt","r")
    bb = open("_q_p_info/_q1_info.txt","r")
    cc = open("_q_p_info/_p0_info.txt","r")
    dd = open("_q_p_info/_p1_info.txt","r")
   
    ndof, ntraj = params["ndof"], params["ntraj"]
    q, p = MATRIX(ndof,ntraj), MATRIX(ndof,ntraj)

    count = 0
    for line in aa:        
        a = line.strip().split(" ")
        q.set(0,count,float(a[0]))
        count += 1

    count = 0
    for line in bb:
        b = line.strip().split(" ")
        q.set(1,count,float(b[0]))
        count += 1

    count = 0
    for line in cc:
        c = line.strip().split(" ")
        p.set(0,count,float(c[0]))
        count += 1

    count = 0
    for line in dd:
        d = line.strip().split(" ")
        p.set(1,count,float(d[0]))
        count += 1

    return q, p

def traj_counter(q, barrier, dof):
    """
    q - MATRIX of trajectory posititions of sinze (ndof x ntraj)
    barrier - speciffied position coordiate. trajectories beyond 
              this point are considered to have "tunneled"
    dof - the degree of freedom to which the coordinate describing
          the barrier is associated with
    Returns the number of trajectories past a speciffied position coordiate,
    for a specifced degree of freedom
    """

    ntraj = q.num_of_cols

    count = 0.0
    for traj in xrange(ntraj):

        if q.get(dof,traj) > 0:          

            count += 1.0


    count /= ntraj
 
    return count

def hermite(n, x):

    r,s,t = 0.0, 0.0, 0.0
    p,q = 1.0, 0.0

    for m in xrange(n):
        r,s = p,q
        p = 2.0*x*r - 2.0*m*t
        q = 2.0*(m+1)*r
        t = r

    return p


def ket_n(q, n, k, m):
    """
    HO state |n>
    """

    hbar = 1.0  # atomic units
    omega = math.sqrt(k/m)  
    alp = m*omega/hbar

    N_n =  math.pow(alp/math.pi, 0.25) / math.sqrt(math.pow(2.0, n) * FACTORIAL(n))
    ksi = math.sqrt(alp)*q    
    H_n = hermite(n, ksi)
 
    res = N_n * H_n * math.exp(-0.5*ksi*ksi)

    return res      


def HO_sup(q, params):
    """
    The probability density function: superposition of HO eigenstates

    """

    k = params["k"]
    m = params["m"]
    states = params["states"]
    coeffs = params["coeffs"]

    x = q.get(0)

    sz = len(states)
    p = 0.0
    for n in xrange(sz):
        p = p + coeffs[n] * ket_n(x, states[n], k, m)

    p = p * p 

    return p

