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



    dt, x, p0, m = 0.1, -2.0, 60.0, 2000.0

    q = MATRIX(1,1); q.set(0, 0, x)
    p = MATRIX(1,1); p.set(0, 0, p0)
    f = MATRIX(1,1); 
    iM = MATRIX(1,1); iM.set(0, 0, 1.0/m)
    Hvib = CMATRIX(2,2)

    Uprev = CMATRIX(2,2)
    Uprev.identity()

    Q = Electronic(2)

    perm_cum = Py2Cpp_int([0, 1])

    while x < 2.5 and x >-2.5:      

        C2El(Cadi, Q)
        propagate_electronic(0.5*dt, Q, Hvib);
        El2C(Q, Cadi)


        p = p + f*dt*0.5 

        q = q + (dt/m)*p
        ham.compute_diabatic(compute_model, q, {})
        ham.compute_adiabatic(1, 0);       

        #######################
        #  Reordering
        X = Uprev * U.H()        
        perm_t = get_reordering(X)
        ham.update_ordering(perm_t)

        perm_cum = composite_permutation(perm_t, perm_cum)
        print Cpp2Py(ham.get_ordering_adi()), Cpp2Py(perm_cum)


        states.permute_rows(perm_t)
        state = perm_t[state]
        Cadi.permute_rows(perm_t)


        """

        u = CMATRIX(U)
        c = CMATRIX(Cadi)
        s = CMATRIX(states)
        hadi = CMATRIX(Hadi)
        fadi = CMATRIX(d1ham_adi[0])
        dcadi = CMATRIX(dc1_adi[0])
        nacadi = CMATRIX(NACadi)
        hvibadi = CMATRIX(Hvib_adi)


        print Cpp2Py(perm_cum)
        print "U_before = "; U.show_matrix()
        print "Hadi_before = "; Hadi.show_matrix()
        push_submatrix(U, u, Py2Cpp_int([0,1]), perm_cum)
        push_submatrix(Hadi, hadi, perm_cum, perm_cum)
        push_submatrix(d1ham_adi[0], fadi, perm_cum, perm_cum)
        push_submatrix(dc1_adi[0], dcadi, perm_cum, perm_cum)
#        push_submatrix(NACadi, nacadi, perm, perm)
#        push_submatrix(Hvib_adi, hvibadi, perm, perm)

#        push_submatrix(Cadi, c, perm_t, Py2Cpp_int([0]))
#        push_submatrix(states, s, perm_t, Py2Cpp_int([0]))
#        state = perm_t[state]

        print "U_after = "; U.show_matrix()
        print "Hadi_after = "; Hadi.show_matrix()

        """
      
#        X = Uprev * U.H()        
        Uprev = CMATRIX(U)
        #######################

        ham.compute_nac_adi(p, iM);
        ham.compute_hvib_adi();       
        Hvib = ham.get_hvib_adi();


        f = ham.forces_adi(states).real().get(0,0)
        p = p + f*dt*0.5 


        ham.compute_nac_adi(p, iM);
        ham.compute_hvib_adi(); 
        Hvib = ham.get_hvib_adi();

        C2El(Cadi, Q)
        propagate_electronic(0.5*dt, Q, Hvib);
        El2C(Q, Cadi)
        

                
        g = compute_hopping_probabilities_fssh(Cadi, ham, 1, dt, 0, 300.0);
        ksi = rnd.uniform(0.0,1.0); 
        new_state = hop(state, g, ksi); 
        new_state = rescale_velocities_adiabatic(p, iM, ham, new_state, state, 1);
        states.set(state, 0, 0.0+0.0j)
        states.set(new_state, 0, 1.0+0.0j)
        state = new_state
 
        dm = Cadi * Cadi.H() 
        x = q.get(0)
        print x, dm.get(0,0).real, dm.get(1,1).real, X.get(0,0).real, X.get(1,1).real, state,\
              Hadi.get(0,0).real, Hadi.get(1,1).real, Hvib.get(0,1).real, Hvib.get(0,1).imag
        
              

#        H = ham.get_ham_dia()
#        S = ham.get_ovlp_dia()
#        check_reordering(H,S, 1.0)
#        Uprev = CMATRIX(U)




run_test()

        
