from liblibra_core import *
from libra_py import *
import math
import sys

def compute_energy(n,omega):

    E = omega * ( n + 0.5 )
    return E


def compute_q(params,coeff_curr,energy,t):

    q = 0.0
    mass  = params["mass"]
    omega = params["omega"]
    prefix = math.sqrt(0.5/(mass*omega))
    nbasis = coeff_curr.num_of_rows

    for n in xrange(nbasis - 1):

        c1 = coeff_curr.get(n+1,0).conjugate() * coeff_curr.get(n,0)
        c2 = coeff_curr.get(n,0).conjugate()   * coeff_curr.get(n+1,0)

        e1 = (energy.get(n,0)-energy.get(n+1,0))
        e2 = (energy.get(n+1,0)-energy.get(n,0))

        exp1 = math.cos( t*e1 ) - 1.0j*math.sin( t*e1 )
        exp2 = math.cos( t*e2 ) - 1.0j*math.sin( t*e2 )

        q +=  math.sqrt(n+1) * ( c1*exp1 + c2*exp2 )

    q *= prefix    

    return q.real


def compute_q2(params,coeff_curr,energy,t):

    q = 0.0
    mass  = params["mass"]
    omega = params["omega"]
    prefix = 0.5/(mass*omega)
    nbasis = coeff_curr.num_of_rows

    # for ( a_+ a_- ) + ( a_- a_+ )
    for n in xrange(nbasis):

        c1 = coeff_curr.get(n,0).conjugate() * coeff_curr.get(n,0)
        q += c1 * (2*n+1)

    # for ( a_+^2 ) + ( a_-^2 )
    for n in xrange(nbasis - 2):

        c2   = coeff_curr.get(n+2,0).conjugate() * coeff_curr.get(n,0)
        c3   = coeff_curr.get(n,0).conjugate() * coeff_curr.get(n+2,0)

        e2   = (energy.get(n,0)-energy.get(n+2,0))  
        e3   = (energy.get(n+2,0)-energy.get(n,0))

        exp2 = math.cos( t*e2 ) - 1.0j*math.sin( t*e2 )
        exp3 = math.cos( t*e3 ) - 1.0j*math.sin( t*e3 )
 
        q += math.sqrt(n+1) * math.sqrt(n+2) * ( c2*exp2 + c3*exp3 )


    q *= prefix
    return q.real



def compute_p(params,coeff_curr,energy,t):

    p = 0.0
    mass  = params["mass"]
    omega = params["omega"]
    prefix = 1.0j*math.sqrt(0.5*mass*omega)
    nbasis = coeff_curr.num_of_rows

    for n in xrange(nbasis - 1):

        c1 = coeff_curr.get(n+1,0).conjugate() * coeff_curr.get(n,0)
        c2 = coeff_curr.get(n,0).conjugate() * coeff_curr.get(n+1,0)

        e1 = (energy.get(n,0)-energy.get(n+1,0))
        e2 = (energy.get(n+1,0)-energy.get(n,0))

        exp1 = math.cos( t*e1 ) - 1.0j*math.sin( t*e1 )
        exp2 = math.cos( t*e2 ) - 1.0j*math.sin( t*e2 )

        p +=  math.sqrt(n+1) * ( c1*exp1 - c2*exp2 )
 

    p *= prefix
    return p.real



def compute_p2(params,coeff_curr,energy,t):

    p = 0.0
    mass  = params["mass"]
    omega = params["omega"]
    prefix = 0.5*(mass*omega)
    nbasis = coeff_curr.num_of_rows

    # for ( a_+ a_- ) + ( a_- a_+ )
    for n in xrange(nbasis):

        c1 = coeff_curr.get(n,0).conjugate() * coeff_curr.get(n,0)
        p +=  c1.real * (2*n+1)

    # for ( a_+^2 ) + ( a_-^2 )
    for n in xrange(nbasis - 2):

        c2   = coeff_curr.get(n+2,0).conjugate() * coeff_curr.get(n,0)
        c3   = coeff_curr.get(n,0).conjugate() * coeff_curr.get(n+2,0)

        e2   = (energy.get(n,0)-energy.get(n+2,0))
        e3   = (energy.get(n+2,0)-energy.get(n,0))

        exp2 = math.cos( t*e2 ) - 1.0j*math.sin( t*e2 )
        exp3 = math.cos( t*e3 ) - 1.0j*math.sin( t*e3 )

        p -= math.sqrt(n+1) * math.sqrt(n+2) * ( c2*exp2 + c3*exp3 )


    p *= prefix
    return p.real


def main(nsnaps, nsteps, params, case):    

    dt = params["dt"]
    coeff = params["coeff"]
    etot = 0.0

    nbasis = len(coeff)
    # Initialize the coeffcient list
    coeff_orig = CMATRIX(nbasis,1)
    coeff_curr = CMATRIX(nbasis,1)
    energy = MATRIX(nbasis,1)
    for n in xrange(nbasis):
        coeff_orig.set(n,0,coeff[n])
        energy.set(n,0,compute_energy(n,params["omega"]))
        etot += compute_energy(n,params["omega"])

    f = open("_res"+str(case)+".txt", "w");  f.close()
    for i in xrange(nsnaps):
        for j in xrange(nsteps):
       
            # For each time step, we have to update all coeffcients 
            for n in xrange(nbasis):

                exp = math.cos( i*(n+0.5) ) - 1.0j*math.sin( i*(n+0.5 ) )
                coeff_curr.set(n,0,coeff_orig.get(n,0)*exp)
            
        # Compute propeties (expectation values, etc)
        q  = compute_q(params,coeff_curr,energy,i)
        p  = compute_p(params,coeff_curr,energy,i)
        q2 = compute_q2(params,coeff_curr,energy,i)
        p2 = compute_p2(params,coeff_curr,energy,i)
   
        # Print results
        f = open("_res"+str(case)+".txt", "a")
        f.write("%8.5f  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" % (i*dt*nsteps, etot, q, p, q2, p2) )
        f.close()



def init_params(case):
    """
    This function intiializes the params dictionary for a given case    
    """

    params = {"mass":2000.0, "omega":0.004, "dt":1.0}

    if   case == 0:
        params.update( { "coeff":[1.0] } )
    elif case == 1:
        params.update( { "coeff":[1.0, 1.0] } )
    elif case == 2:
        params.update( { "coeff":[1.0, 1.0, 1.0] } )

    return params


def run_HO(nsnaps, nsteps):
    for case in [0,1,2]:

        params = init_params(case)
        main(nsnaps, nsteps, params, case)

nsnaps = 50
nsteps = 10
run_HO(nsnaps, nsteps)
