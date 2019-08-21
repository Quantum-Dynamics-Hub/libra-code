import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import units
from libra_py import fit


def lattice_hop(r, T, params, rnd):
    """
    Args:
        r (VECTOR): original coordinate [Bohr]
        T (double): temperature [K]
        params (dictionary): hopping model parameters
       
            * a (VECTOR): lattice vector 1 [Bohr]
            * b (VECTOR): lattice vector 2 [Bohr]
            * c (VECTOR): lattice vector 3 [Bohr]
            * Ea (double): activation energy for site hop [a.u.]

        rnd (Random): random numbers generator object
        
    Returns: 
        VECTOR: new coordinate

    """

    model = params["model"]

    R = None

    if model==1 or model==2:  # regular 2D diffusion

        a,b,c = params["a"], params["b"], params["c"]
        Ea = params["Ea"]

        ksi = rnd.uniform(0.0, 1.0)
        direct = rnd.uniform(0.0, 4.0)   # 2D diffusion
        prob = math.exp(-Ea/(units.kB * T))

        if ksi<prob:
            if direct<=1.0:
                R = r + a
            elif direct<=2.0 and direct > 1.0:
                R = r - a
            elif direct<=3.0 and direct > 2.0:
                R = r + b
            elif direct<=4.0 and direct > 3.0:
                R = r - b
            elif direct<=5.0 and direct > 4.0:
                R = r + c
            elif direct<=6.0 and direct > 5.0:
                R = r - c
        else:
            R = r

    elif model==2:  # additional random displacement

        diam = params["diam"]
        Ea2 = params["Ea2"]    

        ksix = rnd.uniform(0.0, 1.0)
        ksiy = rnd.uniform(0.0, 1.0)
        ksir = math.sqrt(ksix**2 + ksiy**2)
        if ksir>0.0:
            ksix = ksix / ksir
            ksiy = ksiy / ksir

        ksi = rnd.uniform(0.0, 1.0)
        prob1 = math.exp(-Ea2/(units.kB * T))

        if ksi<prob1:
            R.x = R.x + diam*ksix
            R.y = R.y + diam*ksiy

    return R






def run_diffusion(T, rnd, params, nsteps, nat):

    # Statistics
    r = VECTORList()
    for at in range(0, nat):
        r.append(VECTOR(0.0, 0.0, 0.0))
 
    time, MSD, D = [], [], []
    t = 0
    for n in range(0, nsteps):
        # Compute the msd
        msd = 0.0        
        for at in range(0, nat):
            msd = msd + r[at].length2()
        msd = msd/float(nat)
        
        time.append(t)
        MSD.append(msd)
        D.append(msd/(t+1.0))

        # Do all the hops
        for at in range(0, nat):
            r[at] = lattice_hop(r[at], T, params, rnd)

        t = t + 1.0

    return time, MSD, D



def run_T_scan(T_vals, rnd, params, nsteps, nat):

    A, B = [], []
    # Run the diffusion
    for T in T_vals:
        time, MSD, D = run_diffusion(T, rnd, params, nsteps, nat)

        # Process data
        sz = len(time)
        ln_time, ln_MSD = [], []
        for i in range(0, sz):
            if time[i]>0.0 and MSD[i]>0.0:
                ln_time.append(math.log(time[i]))
                ln_MSD.append(math.log(MSD[i]))

        a, b = None, None
        if len(ln_time)>1:
            a, b = fit.Regression(ln_time, ln_MSD, 1)        
        A.append(a)
        B.append(b)

    return A, B
  

