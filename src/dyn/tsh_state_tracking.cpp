/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_state_tracking.cpp
  \brief The file implements the algorithms for state tracking
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "../util/libutil.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libutil;




vector<int> get_reordering(CMATRIX& time_overlap){
    /**
    """ This function identifies which states have changed their identities during some
    calculations (usually the eigenvalue problem) in comparison to what they might have been.

    In the dynamics, this situation occurs when the system passes the conical intersection region
    and the identity of the states may change. The states' energies become non-descriptive for this
    purpose and one needs to look at the changes of the orbitals (e.g. eigenvectors). In particular,
    we can look at the overlap of the sates at adjacent times: <phi_i(t)|phi_i(t+dt)> 
    If no spurious state changes happens, the diagonal elements should be close to 1.0. 
    If they are not - we locate to which state the transitions might have happened.

    In general context, the "time_overlap" matrix is compused of the overlaps of the eigenvectors
    for two problems - the original one and a perturbed one.

    \param[in] time_overlap ( CMATRIX ) the time overlap matrix, <phi_i(t)|phi_j(t+dt)>.

    Returns:
    perm - list of integers that describe the permutation. That is:
    perm[i] - is the index identifying the "older" state "i". Now, it may be labeled
    by some other index, j = perm[i]. 

    """
    */

    // extract the indices where <phi_i(t)|phi_i(t+dt)> is not close to 1. 
    CMATRIX S(time_overlap);  // just a temporary working object
    int sz = time_overlap.n_rows;
    int i;

    // Original permutation
    vector<int> perm(sz, 0);  
    vector<int> perm_cum(sz, 0);  
    for(i=0;i<sz;i++){ perm_cum[i] = i; } 
    
    for(int col=0; col<sz; col++){

      int indx = -1;
      complex<double> val(0.0, 0.0);
      
      while(indx!=col){

        // Find the max element in the given column "col"
        S.max_col_elt(col, val, indx);
            
        // Apply the permutation (col, indx) to the present "perm" list
        for(i=0;i<sz;i++){ perm[i] = i; } 

        int tmp = perm[col];
        perm[col] = perm[indx];
        perm[indx] = tmp;

        // Do the corresponding swap of the columns in the S matrix
        S.swap_cols(col,indx);

        update_permutation(perm, perm_cum);


      }// while indx!=col
    }// for col

    return perm_cum;
}


        
         
MATRIX make_cost_mat(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha){
    /**

    Makes the cost matrix from a given TDM and information on states' energies

    Args:    
        orb_mat_inp ( MATRIX(nstates, nstates) ) : the transition density matrix TDM in a given basis, <i(t)|j(t+dt)> 

        en_mat_inp ( MATRIX(nstates, nstates) ) : the energies of states in a given basis [units: a.u.]

        alpha ( double ): Parameter controlling the range of the orbitals that can participate in
            the reordering. Setting is to 0 makes all orbitals be considered for reordering
            Setting it to a large number makes the effective number of orbitals participating
            in the reordering smaller - this can be used to turn off the reordering. [units: a.u.^-1]

    Returns: 
        MATRIX(nstates, nstates): the matrix of the cost values for different pairs of states

    */

    int a,b;
    int nstates = orb_mat_inp.n_cols;
    MATRIX cost_mat(nstates, nstates);

    for(a=0;a<nstates;a++){
        for(b=0;b<nstates;b++){

            complex<double> s = orb_mat_inp.get(a,b);
            double s2 = (s*std::conj(s)).real();
            double dE = (en_mat_inp.get(a,a) - en_mat_inp.get(b,b)).real();

            double val = s2 * exp(-(alpha*dE)*(alpha*dE));
            cost_mat.set(a,b,val);

        }
    }

    return cost_mat;
}

vector<int> Munkres_Kuhn(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha, int verbosity){

    MATRIX cost_mat = make_cost_mat(orb_mat_inp, en_mat_inp, alpha);

    // Solve the optimal assignment problem for diagonal blocks
    return Munkres_Kuhn_maximize(cost_mat, verbosity);

}




vector<int> get_stochastic_reordering(CMATRIX& time_overlap, Random& rnd){
    /**
    """ This function identifies which states have changed their identities during some
    calculations (usually the eigenvalue problem) in comparison to what they might have been.

    In the dynamics, this situation occurs when the system passes the conical intersection region
    and the identity of the states may change. The states' energies become non-descriptive for this
    purpose and one needs to look at the changes of the orbitals (e.g. eigenvectors). In particular,
    we can look at the overlap of the sates at adjacent times: <phi_i(t)|phi_i(t+dt)> 
    If no spurious state changes happens, the diagonal elements should be close to 1.0. 
    If they are not - we locate to which state the transitions might have happened.

    In general context, the "time_overlap" matrix is compused of the overlaps of the eigenvectors
    for two problems - the original one and a perturbed one.

    \param[in] time_overlap ( CMATRIX ) the time overlap matrix, <phi_i(t)|phi_j(t+dt)>.

    Returns:
    perm - list of integers that describe the permutation. That is:
    perm[i] - is the index identifying the "older" state "i". Now, it may be labeled
    by some other index, j = perm[i]. 

    |phi_j(t+dt)> = a_0j |phi_0(t)> + a_1j |phi_1(t)> + ... a_Nj |phi_N(t)>

    So:  <phi_i(t)|phi_j(t+dt)> = a_ij and then
 
    |a_ij|^2 is the probability that the state i at time t becomes state j at time t+dt

    i - state index at t
    perm[i] - state index at t+dt

    """
    */

    // extract the indices where <phi_i(t)|phi_i(t+dt)> is not close to 1. 
    CMATRIX S(time_overlap);  // just a temporary working object
    int sz = time_overlap.n_rows;
    int i,j, target;

    // Original permutation
    vector<int> perm(sz, 0);  
    for(i=0;i<sz;i++){ perm[i] = i; } 


    vector<int> used;  // already utilized options
    MATRIX g(1, sz);   // switching probabilities

    for(i=0;i<sz-1;i++){ // all starting states (at time t)

        // Make a list of the states that can still participate
        // in the reordering
        vector<int> active_states;
        for(j=0;j<sz;j++){ 
            if(is_in_vector(j,used)){  ;;  }
            else{ active_states.push_back(j); }
        }


        if(sz-i>1){
 
            MATRIX g(1,sz-i);
            for(j=0;j<sz-i;j++){ // all target states

                int indx = active_states[j];
                g.set(0,j, (S.get(i,indx)*std::conj(S.get(i,indx))).real() );
            }
 
            /*
            // Exclude those which we have already assigned
            if(is_in_vector(j,used)){   
                g.set(0,j, 0.0); 
            }

            // For others - compute probabilities. Note: the probabilities will
            // be normalized in the hop function, so at this point we only need the
            // relative magnitudes
            else{  
                g.set(0,j, (S.get(i,j)*std::conj(S.get(i,j))).real() );      
            }
            */

        //}// for j

            // Now stochastically determine the states re-assignment
            double ksi = rnd.uniform(0.0, 1.0);
            target = hop(0, g, ksi);    // internal (active states') state index
            target = active_states[target]; // actual state index
        }
        else if(sz-i==1){
            target = active_states[0];
        }

        // Update the permutaiton and add the found target state to the exclude list
        perm[i] = target;
        used.push_back(target);

    }// for i

    for(i=0;i<sz;i++){ 
      if(!is_in_vector(i,perm)){
        for(j=0;j<sz;j++){ 
            cout<<"j= "<<j<<" perms[j] = "<<perm[j]<<endl; 
            exit(0);
        }
      }
    } 

    return perm;
}


}// namespace libdyn
}// liblibra


