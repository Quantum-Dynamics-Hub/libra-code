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

from numpy import *
from numpy.linalg import *
from cmath import *



"""
  This file demonstrates how to run TSH calculations for a single trajectory

"""

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *




def compute_model(q, params, full_id):

    res = models_Tully.Tully1(q, params)
#    res.rep = params["rep"]    

    return res
    

def signature(X):
    I = CMATRIX(X)
    I *= 0.0
    I.identity()

    res = []   
    tmp = I.H() * X;
    ncols = X.num_of_cols

    for n in xrange(ncols):
        val = tmp.get(n,n).real;
        if val>0.0:
            res.append(1)
        else:
            res.append(-1)

    return res


def check_reordering(H,S, dx):

    n = H.num_of_cols
    E = CMATRIX(n,n)
    U = CMATRIX(n,n)
    Ep = CMATRIX(n,n)
    Up = CMATRIX(n,n)
    Hp = CMATRIX(H)

    Hp.add(0,0, -dx*(1.0+0.0j))
    
    solve_eigen(H,S,E,U,0)

#    if(abs(E.get(0,0).real-E.get(1,1).real)<1e-7):

    solve_eigen(Hp,S,Ep,Up,0)

    Map = U * Up.H()
    print unavoided.get_reordering(Map)
    Map.show_matrix()



def C2El(C, el):
    el.q[0], el.p[0] = C.get(0).real, C.get(0).imag
    el.q[1], el.p[1] = C.get(1).real, C.get(1).imag

def El2C(el, C):
    C.set(0, el.q[0]+1.0j*el.p[0])
    C.set(1, el.q[1]+1.0j*el.p[1])


def composite_permutation(perm_t, perm_cum):
    """
    perm_cum - Cumulative permutation
    perm_t - permutation at a given point

    This function computes a composition of 
    the two permutations: perm_t (x) perm_cum

    E.g. if perm_cum = [1,0], meaning perm_cum(0) = 1, perm_cum(1) = 0
    and perm_t = [0,1], meaning perm_t(0) = 0, perm_t(1) = 1
    Then:
    perm_t (x) perm_cum = [1, 0]
    (perm_t * perm_cum)(0) = perm_t( perm_cum(0) ) = perm_t(1) = 1
    (perm_t * perm_cum)(1) = perm_t( perm_cum(1) ) = perm_t(0) = 0
    """

    tmp = intList()
    sz_cum = len(perm_cum)
    for i in xrange(sz_cum):
        tmp.append(perm_t[ perm_cum[i] ])

    return tmp
        

def phase_changes(X, X_prev):

    sz = X.num_of_cols
    nr = X.num_of_rows
    phases = [1.0, 1.0]
    for c in xrange(sz):
        x = X.col(c)
        x_prev = X_prev.col(c)

        f = (x.H() * x_prev).get(0)
        af = abs(f)
        if af > 0.0:
            phases[c] = f / af

    return phases                 

#        for r in xrange(nr):
#            X.scale(r,c, f.conjugate())
#        C.scale(r, 0, f)
#    print phases

def phase_correct1(u, nac_adi, hvib_adi, phases):

    nst = u.num_of_cols

    for i in xrange(nst):
        for j in xrange(nst):
            u.scale(j, i, phases[i].conjugate())

            fij = phases[j] * phases[i].conjugate()
            nac_adi.scale(j, i, fij)
            hvib_adi.scale(j, i, fij)


def phase_correct2(cadi, phases, phases_prev):

    nst = cadi.num_of_rows

    for i in xrange(nst):
        scl = 1.0+0.0j
        if abs(phases_prev[i])>0.0:
            scl = phases[i]/phases_prev[i]
        cadi.scale(i, 0, scl)

       

def run_test():
    rnd = Random()

    I = CMATRIX(2,2)  # internal reference
    I.set(0,0, 1.0+0.0j);     I.set(0,1, 0.0+0.0j);
    I.set(1,0, 0.0+0.0j);     I.set(1,1, 1.0+0.0j);


    Cadi = CMATRIX(2,1)
    Cadi.set(0, 0, 1.0+0.0j)
    Cadi.set(1, 0, 0.0+0.0j)

    state = 0
    states = CMATRIX(2,1)
    states.set(state, 0, 1.0+0.0j)

    ham = nHamiltonian(2,2,1)
#    ham.init_all(2)

    Hdia = CMATRIX(2,2);     ham.set_ham_dia_by_ref(Hdia);    
    Sdia = CMATRIX(2,2);     ham.set_ovlp_dia_by_ref(Sdia); 
    Hadi = CMATRIX(2,2);     ham.set_ham_adi_by_ref(Hadi);  
    NACadi = CMATRIX(2,2);   ham.set_nac_adi_by_ref(NACadi);
    NACdia = CMATRIX(2,2);   ham.set_nac_dia_by_ref(NACdia);
    Hvib_adi = CMATRIX(2,2); ham.set_hvib_adi_by_ref(Hvib_adi);
    Hvib_dia = CMATRIX(2,2); ham.set_hvib_dia_by_ref(Hvib_dia);
    U = CMATRIX(2,2);        ham.set_basis_transform_by_ref(U); 

    d1ham_dia = CMATRIXList(); d1ham_adi = CMATRIXList()
    dc1_dia = CMATRIXList();   dc1_adi = CMATRIXList()
    for i in [0]:
        d1ham_dia.append( CMATRIX(2,2) ); d1ham_adi.append( CMATRIX(2,2) )
        dc1_dia.append( CMATRIX(2,2) );   dc1_adi.append( CMATRIX(2,2) )
    
    ham.set_d1ham_dia_by_ref(d1ham_dia);  ham.set_d1ham_adi_by_ref(d1ham_adi)
    ham.set_dc1_dia_by_ref(dc1_dia);      ham.set_dc1_adi_by_ref(dc1_adi)



    dt, x, p0, m = 0.1, -2.0, 10.0, 2000.0

    q = MATRIX(1,1); q.set(0, 0, x)
    p = MATRIX(1,1); p.set(0, 0, p0)
    f = MATRIX(1,1); 
    iM = MATRIX(1,1); iM.set(0, 0, 1.0/m)
    Hvib = CMATRIX(2,2)

    Uprev = CMATRIX(2,2)
    Uprev.identity()

    ham.compute_diabatic(compute_model, q, {})
    ham.compute_adiabatic(1, 0);   
    

    Q = Electronic(2)

    cum_phases_prev = [1.0+0.0j, 1.0+0.0j]

    perm_cum = Py2Cpp_int([0, 1])

    while x < 2.5 and x >-2.5:      
#        print "== ==", q.get(0)

        int_state = ham.get_ordering_adi()[state]
#        print state, int_state

        ham.compute_nac_adi(p, iM)
        ham.compute_hvib_adi()



#        propagate_electronic(0.5*dt, Cadi, ham, 1)

        Hvib = ham.get_hvib_adi()
        C2El(Cadi, Q)
        propagate_electronic(0.5*dt, Q, Hvib)
        El2C(Q, Cadi)


        states *= 0.0
        states.set(int_state, 0, 1.0+0.0j)
        p = p + ham.Ehrenfest_forces_adi(states, 0).real() *dt*0.5 
#        p = p + ham.Ehrenfest_forces_adi(Cadi, 0).real() *dt*0.5 
        q = q + (dt/m)*p

        Uprev = CMATRIX(U)
        ham.compute_diabatic(compute_model, q, {})
        ham.compute_adiabatic(1, 0);       

       

        ##============ Get the eigenvalues using Numpy ====================
        h = matrix([ [Hdia.get(0,0), Hdia.get(0,1)], [Hdia.get(1,0), Hdia.get(1,1)] ] )
        Eval,Evec = eig(h)
        _Eval = CMATRIX(2,2)
        _Evec = CMATRIX(2,2)
        Evec = array(Evec)
        
        _Eval.set(0,0, Eval[0]);     _Eval.set(1,1, Eval[1])
        _Evec.set(0,0, Evec[0][0]);  _Evec.set(0,1, Evec[0][1]);
        _Evec.set(1,0, Evec[1][0]);  _Evec.set(1,1, Evec[1][1]);

        #(Hdia * _Evec - _Evec * _Eval).show_matrix()  # this should be zero        
        #print Evec
        #print Eval
        #print h

#        print phases

#        print "Numpy eigenvectors: "; _Evec.show_matrix()
#        print "Eigen eigenvectors: "; U.show_matrix()


        #######################        
        print "Prev U = "; Uprev.show_matrix()
        print "Curr U = "; U.show_matrix()
        #  Reordering
        X = Uprev.H() * U
#        X.show_matrix()
        perm_t = get_reordering(X)
        ham.update_ordering(perm_t)
        Cadi.permute_rows(perm_t)

        print "Reordered U = "; U.show_matrix()

        cum_phases = phase_changes(U, Uprev)
        phase_correct1(U, NACadi, Hvib_adi, cum_phases)
        phase_correct2(Cadi, cum_phases, cum_phases_prev)

        cum_phases_prev = list(cum_phases)
        

        print "Phase corrected U = "; U.show_matrix()


       

        int_state = ham.get_ordering_adi()[state]

        #######################

        states *= 0.0
        states.set(int_state, 0, 1.0+0.0j)
        p = p + ham.Ehrenfest_forces_adi(states, 0).real() *dt*0.5 
#        p = p + ham.Ehrenfest_forces_adi(Cadi, 0).real() *dt*0.5 

        ham.compute_nac_adi(p, iM);
        ham.compute_hvib_adi(); 

#        propagate_electronic(0.5*dt, Cadi, ham, 1)

        Hvib = ham.get_hvib_adi()
        C2El(Cadi, Q)
        propagate_electronic(0.5*dt, Q, Hvib)
        El2C(Q, Cadi)

        

        
        g = compute_hopping_probabilities_fssh(Cadi, ham, 1, dt, 0, 300.0);
        ksi = rnd.uniform(0.0,1.0); 

        new_state = hop(int_state, g, ksi); 
        new_state = rescale_velocities_adiabatic(p, iM, ham, new_state, int_state, 1);

        int_state = new_state   # internal state

        # Convert to the physical state:
        state = inverse_permutation(ham.get_ordering_adi())[ int_state ];
        


                  
        dm = Cadi * Cadi.H() 
        x = q.get(0)
        P = ham.get_ordering_adi()

#        print x, dm.get(0,0).real, dm.get(1,1).real, X.get(0,0).real, X.get(1,1).real, state,\
#              Hadi.get(0,0).real, Hadi.get(1,1).real, Hvib.get(0,1).real, Hvib.get(0,1).imag

        print x, dm.get(P[0],P[0]).real, dm.get(P[1],P[1]).real, X.get(0,0).real, X.get(1,1).real, state,\
              Hadi.get(P[0],P[0]).real, Hadi.get(P[1],P[1]).real, Hvib.get(P[0],P[1]).real, Hvib.get(P[0],P[1]).imag, cum_phases[0], cum_phases[1]
        
              

run_test()

        
