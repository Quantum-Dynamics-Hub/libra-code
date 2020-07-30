/*********************************************************************************
* Copyright (C) 2018-2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_methods_dish.cpp
  \brief The file implements the decoherence-induced surface hopping (DISH) method

*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"
#include "dyn_decoherence.h"


/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


vector<int> decoherence_event(MATRIX& coherence_time, MATRIX& coherence_interval, Random& rnd){
/**
  For each trajectory, check which adiabatic states have evolved longer than decoherence times.
  In the case that multiple states have done this, we select only one randomly.

  coherence_time - MATRIX(nst, ntraj) - for each state is how long has that state resided in a coherence evolution
  coherence_interval - MATRIX(1, ntraj) - for each trajectory - what is the longest coherence time in the system

  Return:
  For each trajectory:
  The selection of the states which will experience decoherence events (collapse or projection out)
  -1 means that no basis states have experienced the decoherence event for this trajectory

*/

  int i,traj;
  int ntraj = coherence_time.n_cols;
  int nst = coherence_time.n_rows;

  /// By default, set the indices of the states that "experience" decoherence events to -1
  /// In the analysis, if we encounter the index of -1, we'll know that no decoherence has happened
  /// and will simply continue the coherent evolution.
  vector<int> res(ntraj, -1);

  for(traj=0; traj < ntraj; traj++){

    vector<int> which_decohere; /// which basis states shall experience the decoherence event

    /// Determine all basis states that may decohere at this point
    for(i=0;i<nst;i++){

      /// The state i has evolved coherently for longer than the coherence interval
      /// so it has to experience a decoherence event
      if(coherence_time.get(i, traj) >= coherence_interval.get(i, traj) ) {
        which_decohere.push_back(i);
      }
    }

    if(which_decohere.size()>0){

      /// Select only one of the basis states that will be experience the decoherence event.
      vector<int> selection(1,0);
      randperm(which_decohere.size(), 1, selection);
      res[traj] = selection[0];
    }

  }// for traj

  return res;

}


vector<int> dish_hop_proposal(vector<int>& act_states, CMATRIX& Coeff,
  MATRIX& coherence_time, vector<MATRIX>& decoherence_rates, Random& rnd){

    int i,traj;
    int nst = Coeff.n_rows;
    int ntraj = Coeff.n_cols;


    /// Update coherence intervals
    MATRIX coherence_interval(nst, ntraj); // for DISH
    coherence_interval = coherence_intervals(Coeff, decoherence_rates);

    /// Determine which states may experience decoherence event
    /// If the decohered_states[traj] == -1, this means no basis states on the trajectory traj
    /// have experienced the decoherence event, othervise the variable will contain an index of
    /// such state for each trajectory
    vector<int> decohered_states( decoherence_event(coherence_time, coherence_interval, rnd) );


    /// By default, the proposed states are assumed to be the current ones
    vector<int> proposed_states(act_states);


    for(traj=0; traj < ntraj; traj++){

      int istate = decohered_states[traj];

      /// Exclude the situation when no decoherence event have occured (-1)
      /// in those cases we of course do not want to preject out those states
      if(istate>-1){

        /// No matter what happens with the wavefunctions, the cohrence interval
        /// is reset for the decohered state, since the state has experienced decoherence event
        coherence_time.set(decohered_states[traj], traj, 0.0);


        /// Propose new discrete states: if the state that experiences decoherence is selected
        /// with the probability given by its SE population, it is considered a proposed hopping state
        /// If it is not selected - the current active states become the proposed states (hop to the same state)
        /// If no decoherence event has happened for a trajectory, the proposed states are set to be the current ones
        /// (also hop to the same state)
        /// The hops to the same states are generally excluded in the following states - so it would be the
        /// standard coherent evolution without hopping
        double prob = (std::conj(Coeff.get(istate, traj)) * Coeff.get(istate, traj) ).real();
        double ksi = rnd.uniform(0.0, 1.0);

        if(ksi<=prob){
          proposed_states[traj] = decohered_states[traj];
        }
        else{
            /// Project out the decohered states if they aren't selected
            if(decohered_states[traj]!= act_states[traj]){
                project_out(Coeff, traj, decohered_states[traj]);
            }
        }

      }// istate > -1

    }// for traj

    return proposed_states;

}

void dish_project_out_collapse(vector<int>& old_states, vector<int>& proposed_states, vector<int>& new_states,
  CMATRIX& Coeff, MATRIX& coherence_time, int collapse_option){
/**
  To handle projections after the attempted hop:
*/

  int ntraj = old_states.size();

  for(int traj=0; traj < ntraj; traj++){

    if(proposed_states[traj] != old_states[traj]){
      /// Attempted non-trivial hops

      if(new_states[traj] == proposed_states[traj]){
        /// Successfull hop - collapse onto this new state
        collapse(Coeff, traj, new_states[traj], collapse_option);
      }
      else{
        /// Attempted hop was unsuccessful - project out the proposed state
        project_out(Coeff, traj, proposed_states[traj]);
      }

    }
    else{  /// proposed_states[traj] == old_states[traj]
      /// Trivial hop = the proposed transition is to the starting state
      /// Just collapse onto that state

      collapse(Coeff, traj, old_states[traj], collapse_option);

    }

/*
    if(new_states[traj] == old_states[traj]){
        /// No transition occured

      if(proposed_states[traj] == old_states[traj]){
        /// No attempted transition
        ; ;
      }
      else if(proposed_states[traj] != old_states[traj]){
        /// Attempted the hop, but didn't accept it - project out the proposed_state
        project_out(Coeff, traj, proposed_states[traj]);
        coherence_time.set(proposed_states[traj], traj, 0.0);
      }

    }

    else if(new_states[traj] != old_states[traj]){
      /// The transition to a new state has been accepted - collapse wfc
      /// onto this new state
      collapse(Coeff, traj, new_states[traj], collapse_option);
      coherence_time.set(proposed_states[traj], traj, 0.0);
    }
*/

  }// for traj

}


}// namespace libdyn
}// liblibra
