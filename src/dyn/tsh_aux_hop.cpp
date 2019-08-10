/*********************************************************************************
* Copyright (C) 2015-2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file tsh_aux_hop.cpp
  \brief The file implements the hopping procedure needed for the TSH algorithms.
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;


/// libdyn namespace
namespace libdyn{


vector<int> tsh_vec2indx(CMATRIX& states){
/**
  Convert a set of vectors <states> into the corresponding indices  
  res[i] - is the index of the row which contains maximal value in the i-th column
*/

  vector<int> res(states.n_cols);

  int istate;
  complex<double> tmp;

  for(int i=0; i<states.n_cols; i++){
    states.max_col_elt(i, tmp, istate);
    res[i] = istate;
  }

  return res;

}

void tsh_indx2vec(nHamiltonian& ham, CMATRIX& states, vector<int>& res){
/**
  Convert the active state index (physical states) to the vector with occupation numbers
  The reordering accumulated up until now is taken into accout. This is done for
  1 or many trajectories
*/
  states *= 0.0;

  vector<int> tx; tx = id_permutation(states.n_rows);
  vector<int> ty(1,0); 

  if(states.n_cols!=ham.children.size()){
    cout<<"ERROR in void tsh_indx2vec(nHamiltonian& ham, CMATRIX& states, vector<int>& res)\
          The size of the input matrix ("<<states.n_cols<<") and the number of children \
          Hamiltonians ("<<ham.children.size()<<") are not equal\n";
    exit(0);
  }
  if(states.n_cols!=res.size()){
    cout<<"ERROR in void tsh_indx2vec(nHamiltonian& ham, CMATRIX& states, vector<int>& res)\
          The size of the input matrix ("<<states.n_cols<<") and the number of items in \
          the resul permutation ("<<res.size()<<") are not equal\n";
    exit(0);
  }
  

  for(int i=0; i<states.n_cols; i++){

    // For a given trajectory
    CMATRIX state(states.n_rows,1);         

    // Set the active state to where it would have been without the reordering
    state.set(res[i], 0, complex<double>(1.0, 0.0));

    // Apply the cumulative reordering 
    vector<int> M = ham.children[i]->get_ordering_adi();
    state.permute_rows(M);

    // Insert the vector to the corresponding column
    ty[0] = i;
    push_submatrix(states, state, tx, ty);
    
  }

}


void tsh_physical2internal(nHamiltonian& ham, vector<int>& internal, vector<int>& physical){
/**
  Map the indices of the active states from the physical notation to 
  the internal notation (eigenvector ordering), using the permutations accumulated so far.
*/

  int ntraj = ham.children.size();

  for(int i=0; i<ntraj; i++){
 
    vector<int> M = ham.children[i]->get_ordering_adi();
    internal[i] = M[ physical[i] ]; 

//    cout<<"i = "<<i<<" physical[i]= "<<physical[i]<<" internal[i]= "<<internal[i]<<endl;
  }

}

void tsh_internal2physical(nHamiltonian& ham, vector<int>& internal, vector<int>& physical){
/**
  Map the indices of the active states from the internal notation (eigenvector ordering) to the
  physical, using the permutations accumulated so far.
*/

  int ntraj = ham.children.size();

  for(int i=0; i<ntraj; i++){
 
    vector<int> M = ham.children[i]->get_ordering_adi();
    vector<int> iM = inverse_permutation( M );

    physical[i] = iM[ internal[i] ]; 
  }

}

/**
CMATRIX compute_phases(CMATRIX& U, CMATRIX& U_prev){


  X - the |psi(t')>
  Xprev = |psi(t)>,  t' > t



  int i;
  complex<double> f;
  int nc = U.n_cols;

  CMATRIX phases(nc, 1);

  // Default values
  for(i=0;i<nc; i++){  phases.set(i, 0, 1.0, 0.0); }

  // Compute phase corrections  
  for(i=0; i<nc; i++){

    f = (U_prev.col(i).H() * U.col(i) ).get(0);
    double af = abs(f);

    if(af > 0.0){   phases.set(i, 0, f / af);     }

  }// for i

  return phases;

}
*/


void phase_correct_ampl(CMATRIX& C, CMATRIX& cum_phases, CMATRIX& cum_phases_prev){
/** 
  C - the amplitude in the superposition: |Psi> = |Phi_adi> * C_adi
  cum_phases - the phase corrections at this timestep
  cum_phases - the phase corrections at a previous timestep
*/

  int nst = C.n_rows;

  for(int i=0;i<nst; i++){
    complex<double> scl(1.0, 0.0);

    if(abs(cum_phases_prev.get(i,0) ) > 0.0 ){
      scl = cum_phases.get(i,0)/cum_phases_prev.get(i,0);
      C.scale(i, 0, scl);
    }
  }// for i

}

void phase_correct_ampl(CMATRIX* C, CMATRIX* phases){
/** 
  C - the amplitude in the superposition: |Psi> = |Phi_adi> * C_adi
  phases - the phase corrections at this timestep

*/

  int nst = C->n_rows;

  for(int i=0;i<nst; i++){  C->scale(i, 0, phases->get(i,0));  }

}

void phase_correct_ampl(CMATRIX& C, CMATRIX& phases){

  phase_correct_ampl(&C, &phases);

}






int hop(int initstate, MATRIX& g, double ksi){
/** 
  \brief Attempts a stochastic hop from the initial state "initstate"
  \param[in] initstate The index of the state from which we try to hop out 
  \param[in] g The hopping probabilities matrix (type MATRIX)
  \param[in] ksi A random number that determines the outcome of the "hop" procedure

  Returned value: the index of the state to which we have hopped
*/
  int i;
  int nstates = g.n_cols;
  double left, right; left = right = 0.0;
  int finstate = -1;
  

  // To avoid problems, lets renormalize the hopping probabilities
  double nrm = 0.0;
  for(i=0;i<nstates;i++){  nrm += g.get(initstate, i);  }

  if(nrm>0.0){

      for(i=0;i<nstates;i++){    
          if(i==0){left = 0.0; right = g.get(initstate,i)/nrm; }
          else{  left = right; right = right + g.get(initstate,i)/nrm; }

         if((left<=ksi) && (ksi<=right)){  finstate = i;  }    
      }

  }
  else{  finstate = initstate; }  // probability to hop to any other states is zero
                                  // so stay on the original state

  if(finstate==-1){
    std::cout<<"Something is wrong in the hop(...) function\nExiting now...\n";
    exit(0);
  }

  return finstate;

}// hop



void hop(int& initstate, Nuclear* mol, Hamiltonian* ham, double ksi, MATRIX* g, int do_rescaling, int rep, int do_reverse){
/** 
  \brief Do actual hop from the state initstate 
  \param[in,out] initstate The state from which we try to hop out - it will also be updated after the hop has happened
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] ksi A random number that determines the outcome of the "hop" procedure
  \param[in] g The hopping probabilities matrix
  \param[in] do_rescaling The flag to turn on/off CPA: 0 - no velocity rescaling (CPA, no back-reaction),
  in this case one should use Boltzmann factor (consider use_boltz_factor when computing the hopping probability matrix, g)
  1 - do rescaling (back-reaction), in this case it would be wrong to use Boltzmann factor
  \param[in] rep Selects the used representation:  0 - for diabatic, 1 - for adiabatic
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.
*/

  int nstates = g->n_cols;
  double left, right; left = right = 0.0;
  int finstate = initstate;

  for(int i=0;i<nstates;i++){
    if(i==0){left = 0.0; right = g->get(initstate,i); }
    else{  left = right; right = right + g->get(initstate,i); }
 
    if((left<ksi) && (ksi<=right)){  finstate = i;  }
  }


  if(finstate!=initstate){

    if(!do_rescaling){ initstate = finstate; }        // CPA-style, no velocity rescaling
    else{                                             // Possibly rescale velocities - normal inclusion of back-electron reaction

      // state is changed or preserved in the function
      if(rep==0){
        rescale_velocities_diabatic(mol,ham,finstate,initstate); 
      }
      else if(rep==1){
        rescale_velocities_adiabatic(mol,ham,finstate,initstate,do_reverse); 
      }

    }// else

  }// finstate!=initstate

  initstate = finstate;

}// hop

int hop(int initstate, Nuclear& mol, Hamiltonian& ham, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse){
/** 
  \brief Do actual hop from the state initstate - Python-friendly
  \param[in] initstate The state from which we try to hop out - it will also be updated after the hop has happened
  \param[in,out] mol Nuclear DOF. Can be updated (velocity rescaling or reversal)
  \param[in,out] ham A handler of Hamiltonian. Internal parameters may be updated, if the Hamiltonian is recomputed
  \param[in] ksi A random number that determines the outcome of the "hop" procedure
  \param[in] g The hopping probabilities matrix
  \param[in] do_rescaling The flag to turn on/off CPA: 0 - no velocity rescaling (CPA, no back-reaction),
  in this case one should use Boltzmann factor (consider use_boltz_factor when computing the hopping probability matrix, g)
  1 - do rescaling (back-reaction), in this case it would be wrong to use Boltzmann factor
  \param[in] rep Selects the used representation:  0 - for diabatic, 1 - for adiabatic
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  The function returns the index of the final state (new or old).
*/


  int res = initstate; 
  hop(res, &mol, &ham, ksi, &g, do_rescaling, rep, do_reverse);

  return res;

}

int hop(int initstate, Ensemble& ens, int i, double ksi, MATRIX& g, int do_rescaling, int rep, int do_reverse){
/** 
  \brief Do actual hop from the state initstate of the i-th trajectory of an ensemble - Python-friendly
  \param[in] initstate The state from which we try to hop out - it will also be updated after the hop has happened
  \param[in] i The index of the trajectory of interest
  \param[in,out] ens Describes the ensemble of trajectories which we propagate (including hops)
  \param[in] ksi A random number that determines the outcome of the "hop" procedure
  \param[in] g The hopping probabilities matrix
  \param[in] do_rescaling The flag to turn on/off CPA: 0 - no velocity rescaling (CPA, no back-reaction),
  in this case one should use Boltzmann factor (consider use_boltz_factor when computing the hopping probability matrix, g)
  1 - do rescaling (back-reaction), in this case it would be wrong to use Boltzmann factor
  \param[in] rep Selects the used representation:  0 - for diabatic, 1 - for adiabatic
  \param[in] do_reverse The option that determines what to do if the hop was rejected because of the energy conservation
  (frustrated hop): do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta(velocities)
  are inverted.

  The function returns the index of the final state (new or old).
*/


  int res = initstate; 
  hop(res, &ens.mol[i], ens.ham[i], ksi, &g, do_rescaling, rep, do_reverse);

  return res;

}





void hop(int ntraj, vector<int>& initstate, vector<Nuclear*>& mol, vector<Hamiltonian*>& ham, 
         vector<double> ksi, vector<MATRIX*>& g, int do_rescaling, int rep, int do_reverse){
// Do actual hop from state initstate
// initstate - state from which we try to hop out - it will also be updated after the hop has happened
// mol - nuclear DOF
// ham - handler of Hamiltonian
// ksi - a random number
// g   - hopping probabilities matrix
// do_rescaling - flag to turn on/off CPA: 0 - no rescaling (CPA), 1 - do rescaling (back-reaction)
// rep - representation:  0 - for diabatic, 1 - for adiabatic
  int traj;
  int nstates = g[0]->n_cols;
  double left, right; 
  vector<int> finstate; finstate = initstate;

  for(traj=0;traj<ntraj;traj++){

    left = right = 0.0;

    for(int i=0;i<nstates;i++){
      if(i==0){left = 0.0; right = g[traj]->get(initstate[traj],i); }
      else{  left = right; right = right + g[traj]->get(initstate[traj],i); }
 
      if((left<ksi[traj]) && (ksi[traj]<=right)){  finstate[traj] = i;  }

    }// for i

  }// for all copies (trajectories) of the system

  int status = 1; // initial and final are equal
  for(traj=0; traj<ntraj; traj++){  if(finstate[traj]!=initstate[traj]){ status = 0; } }

  if(status==0){  // different multi-particle states

    if(!do_rescaling){ initstate = finstate; }        // CPA-style, no velocity rescaling
    else{                                             // Possibly rescale velocities - normal inclusion of back-electron reaction

      // state is changed or preserved in the function
      if(rep==0){
        //rescale_velocities_diabatic(mol,ham,finstate,initstate); 
      }
      else if(rep==1){
        rescale_velocities_adiabatic(ntraj,mol,ham,finstate,initstate,do_reverse); 
      }

    }// else

  }// finstate!=initstate

  initstate = finstate;

}// hop


vector<int>
hop(int ntraj, vector<int> initstate, vector<Nuclear>& mol, vector<Hamiltonian>& ham, 
    vector<double> ksi, vector<MATRIX>& g, int do_rescaling, int rep, int do_reverse){

/*
  vector<int> res; res = initstate;
  vector<Nuclear*> _mol;
  vector<Hamiltonian*> _ham;
  vector<MATRIX*> _g;

  for(int traj=0; traj<ntraj; traj++){
    _mol[traj] = &mol[traj];
    _ham[traj] = &ham[traj];
    _g[traj] = &g[traj];
  }


  hop(ntraj, res, _mol, _ham, ksi, _g, do_rescaling, rep, do_reverse);

  return res;
*/
  int traj;
  int nstates = g[0].n_cols;
  double left, right; 
  vector<int> finstate; finstate = initstate;

  for(traj=0;traj<ntraj;traj++){

    left = right = 0.0;

    for(int i=0;i<nstates;i++){
      if(i==0){left = 0.0; right = g[traj].get(initstate[traj],i); }
      else{  left = right; right = right + g[traj].get(initstate[traj],i); }
 
      if((left<ksi[traj]) && (ksi[traj]<=right)){  finstate[traj] = i;  }

    }// for i

  }// for all copies (trajectories) of the system

  int status = 1; // initial and final are equal
  for(traj=0; traj<ntraj; traj++){  if(finstate[traj]!=initstate[traj]){ status = 0; } }

  if(status==0){  // different multi-particle states

    if(!do_rescaling){ initstate = finstate; }        // CPA-style, no velocity rescaling
    else{                                             // Possibly rescale velocities - normal inclusion of back-electron reaction

      // state is changed or preserved in the function
      if(rep==0){
        //rescale_velocities_diabatic(mol,ham,finstate,initstate); 
      }
      else if(rep==1){
        rescale_velocities_adiabatic(ntraj,mol,ham,finstate,initstate,do_reverse); 
      }

    }// else

  }// finstate!=initstate


  return finstate;  

}


boost::python::list
hop(int ntraj, boost::python::list initstate, boost::python::list mol, boost::python::list ham, 
    boost::python::list ksi, boost::python::list g, int do_rescaling, int rep, int do_reverse){
/*
  vector<int> res; 
  vector<Nuclear*> _mol;
  vector<Hamiltonian*> _ham;
  vector<MATRIX*> _g;

  for(int traj=0; traj<ntraj; traj++){
     res[traj] = extract<int>(initstate[traj]);
//    _mol[traj] = extract<Nuclear>(mol[traj]);
//    _ham[traj] = extract<Hamiltonian>(ham[traj]);
//    _g[traj] = extract<MATRIX>(g[traj]);
  }

//  hop(ntraj, res, _mol, _ham, ksi, _g, do_rescaling, rep, do_reverse);

//  boost::python::list lres; lres = [];

//  for(int traj=0; traj<ntraj; traj++){
//     lres.append(res[traj]);
//  }

  return lres;
*/   
}




}// namespace libdyn
}// liblibra

