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

void tsh_indx2vec(CMATRIX& states, vector<int>& res){
/**
  Sets occupies the trajectories according to the scheme res
*/
  states *= 0.0;

  for(int i=0; i<states.n_cols; i++){
    states.set(res[i], i, complex<double>(1.0, 0.0));
  }

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

  for(i=0;i<nstates;i++){
    if(i==0){left = 0.0; right = g.get(initstate,i)/nrm; }
    else{  left = right; right = right + g.get(initstate,i)/nrm; }
 
    if((left<ksi) && (ksi<=right)){  finstate = i;  }
  }

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

