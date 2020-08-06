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
  \file dyn_projectors.cpp
  \brief The file implements the algorithms for updating the state projectors
  that account for state tracking and phase corrections


  CMATRIX compute_phase_corrections(CMATRIX& S, double tol)
  CMATRIX compute_phase_corrections(CMATRIX& S)
  vector<int> get_reordering(CMATRIX& time_overlap)
  MATRIX make_cost_mat(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha)
  vector<int> Munkres_Kuhn(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha, int verbosity)
  vector<int> get_stochastic_reordering(CMATRIX& time_overlap, Random& rnd)
  CMATRIX permute2cmatrix(vector<int>& permutation)
  void update_projectors(dyn_control_params& prms, vector<CMATRIX>& projectors,
    vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd)

*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "Dynamics.h"
#include "../util/libutil.h"
#include "dyn_control_params.h"


/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

using namespace libutil;

namespace bp = boost::python;


CMATRIX compute_phase_corrections(CMATRIX& S, double tol){
/**
  Compute the phase correction of one set of eigenvectors (new) with respect to
  another one (reference).

  S = <psi(ref)|psi(new)> - N x N matrix, where N is the number of states

  Result - N x 1 matrix of the phase correction factors.

  To obtain the phase-corrected wavefunction, multiply the new wavefunctions
  by the complex-conjugated factors:

  |psi_i(new, phase-corrected)> = |psi_i(new)> * phase_corr.get(i,0).conjugate()

*/

  int i;
  complex<double> f;
  int nc = S.n_cols;

  CMATRIX phase_corr(nc, 1);

  // Default values
  for(i=0;i<nc; i++){  phase_corr.set(i, 0, 1.0, 0.0); }

  // Compute phase corrections
  for(i=0; i<nc; i++){

    f = S.get(i,i);

    double af = sqrt( std::norm(f) ); // std::norm(val) returns the norm squared of complex number val. So we use sqrt()
    //double af = std::norm(f);


    // If the overlap is very small, it means we might have switched
    // the state adiabatically, in which case the phase correction does
    // not make sense
    if(af > tol){   phase_corr.set(i, 0, f / af);     }

  }// for i

  return phase_corr;

}

CMATRIX compute_phase_corrections(CMATRIX& S){

  return compute_phase_corrections(S, 1e-3);
}





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

    // was sz - 1 - not sure why
    for(i=0;i<sz;i++){ // all starting states (at time t)

        // Make a list of the states that can still participate
        // in the reordering
        vector<int> active_states;
        for(j=0;j<sz;j++){
            if(is_in_vector(j,used)){  ;;  }
            else{ active_states.push_back(j); }
        }


        int act_st_sz = active_states.size(); // should be equal to sz - i

        if(act_st_sz > 1){

            MATRIX g(1, act_st_sz);
            for(j=0; j<act_st_sz; j++){ // all target states

                int indx = active_states[j];
                g.set(0,j, (S.get(i,indx)*std::conj(S.get(i,indx))).real() );
            }

            // Now stochastically determine the states re-assignment
            double ksi = rnd.uniform(0.0, 1.0);

            target = hop(0, g, ksi);    // internal (active states') state index
            if(target > act_st_sz-1){
                cout<<"ERROR: target = "<<target<<" is larger than max possible value = "<<act_st_sz-1<<endl;
                exit(0);
            }
            target = active_states[target]; // actual state index
        }
        else if(act_st_sz == 1){
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

vector<int> get_stochastic_reordering2(CMATRIX& time_overlap, Random& rnd){
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

    int initial_state = 0;
    double ksi = rnd.uniform(0.0, 1.0);
    int number_states = time_overlap.n_rows;
    vector<int> list_states;
    vector<int> permutation;
    vector<double> all_probablities;
    double permutation_probability;
    double sum_overlap;
    double state_probability;

// Construct perumtations and calculate probability of each permutation

    for(int i = 0; i<number_states; i++)
        list_states.push_back(i);

    vector< vector<int> > list_permutations = compute_all_permutations(list_states);
    MATRIX g(1,list_permutations.size());

    for(int permutation_num = 0; permutation_num<list_permutations.size(); permutation_num++){
        permutation_probability = 1;
        sum_overlap = 0;
        permutation = list_permutations[permutation_num];
        for(int i = 0; i<number_states; i++){
                state_probability = real(conj(S.get(i, permutation[i])) * S.get(i, permutation[i]));
                    permutation_probability *= state_probability;
                all_probablities.push_back(permutation_probability);
                sum_overlap += state_probability;
        } // for i
        g.set(0, permutation_num, (permutation_probability/sum_overlap));
    } // for permutation

// normalize probabilities
double sum_probabilities = 0;
for(int probability_index = 0; probability_index<list_permutations.size(); probability_index++)
    sum_probabilities += g.get(0,probability_index);
for(int probability_index = 0; probability_index<list_permutations.size(); probability_index++)
    g.scale(0,probability_index, (1/sum_probabilities));
// Stochastically determine permutations

    vector<int> chosen_permutation = list_permutations[hop(initial_state, g, ksi)];

// just for testing
/*
    double total_prob = 0;
    for(int prob_index = 0; prob_index<list_permutations.size(); prob_index++)
    total_prob += abs(g.get(0,prob_index));


    for(int i=0; i<list_permutations.size(); i++){
        permutation = list_permutations[i];
        cout<< "Probability of permutation ";
        for(int j=0; j<permutation.size(); j++)
            cout << permutation[j];
        cout<< " is "<< g.get(0,i) << endl;
    }
    cout<< "Sum of probabilities is "<< total_prob<< endl<< "The chosen permutation is ";
    for(int k=0; k<permutation.size(); k++)
        cout << chosen_permutation[k];
    cout<< endl;
*/

    return chosen_permutation;

}




CMATRIX permutation2cmatrix(vector<int>& permutation){
/**
  Represent a permutation as a matrix
*/

  int nst = permutation.size();

  CMATRIX res(nst, nst);

  for(int i=0; i<nst; i++){

    res.set(permutation[i], i, complex<double>(1.0, 0.0));

  }

  return res;

}


void update_projectors(dyn_control_params& prms, vector<CMATRIX>& projectors,
  vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd){


  int nst = projectors[0].n_rows;
  int ntraj = projectors.size();

  vector<int> perm_t(nst,0);
  for(int i=0; i<nst; i++){ perm_t[i] = i; }

  CMATRIX phase_i(nst, 1);
  CMATRIX st(nst, nst);
  CMATRIX ist(nst, nst);
  CMATRIX projector_old(nst, nst);

  for(int traj=0; traj<ntraj; traj++){

    projector_old = projectors[traj];
    st = St[traj];

    if(prms.state_tracking_algo==1){
        perm_t = get_reordering(st);
    }
    else if(prms.state_tracking_algo==2){
        perm_t = Munkres_Kuhn(st, Eadi[traj], prms.MK_alpha, prms.MK_verbosity);
    }
    if(prms.state_tracking_algo==3){
        perm_t = get_stochastic_reordering(st, rnd);
    }

    if(prms.state_tracking_algo==32){
        perm_t = get_stochastic_reordering2(st, rnd);
    }

    // P -> P * perm
    CMATRIX p_i(nst, nst);
    p_i = permutation2cmatrix(perm_t);
    projectors[traj] = projectors[traj] * p_i;
    st = projector_old.H() * st * projectors[traj];


    if(prms.do_phase_correction){

      // ### Compute the instantaneous phase correction factors ###
      phase_i = compute_phase_corrections(st);  // f(i)

      // ### Scale projections' components by the phases ###
      for(int a=0; a<nst; a++){
        projectors[traj].scale(-1, a, std::conj(phase_i.get(a)) );
      }
    }


  }// for traj

}




CMATRIX raw_to_dynconsyst(CMATRIX& amplitudes, vector<CMATRIX>& projectors){
/**
  This function converts the amplitudes from the raw ones to the dynamically-consistent

  |dyn_cons> = |raw> * P, so

  |Psi> = |raw> * C_raw = |dyn_cons> * C_dyn_cons

  Then:
  |raw> * C_raw = |raw> * P * C_dyn_cons

   So: C_raw = P * C_dyn_cons
   and
       C_dyn_cons = P.H() * C_raw


*/
  int nst = amplitudes.n_rows;
  int ntraj = amplitudes.n_cols;

  CMATRIX tmp(nst, 1);
  CMATRIX res(nst, ntraj);

  vector<int> x_stenc(1, 0);
  vector<int> y_stenc(nst,0); for(int i=0; i<nst; i++){ y_stenc[i] = i; }

  if(projectors.size() != ntraj){
    cout<<"ERROR in raw_to_dynconsyst: the dimensions do not agree\n";
    exit(0);
  }


  for(int traj = 0; traj < ntraj; traj++){

    x_stenc[0] = traj;
    pop_submatrix(amplitudes, tmp, y_stenc, x_stenc);

    tmp = projectors[traj].H() * tmp;

    push_submatrix(res, tmp, y_stenc, x_stenc);

  }

  return res;

}

CMATRIX dynconsyst_to_raw(CMATRIX& amplitudes, vector<CMATRIX>& projectors){
/**
  This function converts the amplitudes from the dynamically-consistent to raw form

  |dyn_cons> = |raw> * P, so

  |Psi> = |raw> * C_raw = |dyn_cons> * C_dyn_cons

  Then:
  |raw> * C_raw = |raw> * P * C_dyn_cons

   So: C_raw = P * C_dyn_cons
   and
       C_dyn_cons = P.H() * C_raw

*/
  int nst = amplitudes.n_rows;
  int ntraj = amplitudes.n_cols;

  CMATRIX tmp(nst, 1);
  CMATRIX res(nst, ntraj);

  vector<int> x_stenc(1, 0);
  vector<int> y_stenc(nst,0); for(int i=0; i<nst; i++){ y_stenc[i] = i; }

  if(projectors.size() != ntraj){
    cout<<"ERROR in raw_to_dynconsyst: the dimensions do not agree\n";
    exit(0);
  }


  for(int traj = 0; traj < ntraj; traj++){

    x_stenc[0] = traj;
    pop_submatrix(amplitudes, tmp, y_stenc, x_stenc);

    tmp = projectors[traj] * tmp;

    push_submatrix(res, tmp, y_stenc, x_stenc);

  }

  return res;

}





}// namespace libdyn
}// liblibra
