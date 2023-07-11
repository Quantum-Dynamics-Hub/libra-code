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
#include "dyn_hop_acceptance.h"
#include "dyn_hop_proposal.h"
#include "../math_specialfunctions/libspecialfunctions.h"


/// liblibra namespace
namespace liblibra{

using namespace libspecialfunctions;

/// libdyn namespace
namespace libdyn{


vector<int> decoherence_event(MATRIX& coherence_time, MATRIX& coherence_interval, int decoherence_event_option, Random& rnd){
/**
  For each trajectory, check which adiabatic states have evolved longer than decoherence times.
  In the case that multiple states have done this, we select only one randomly.
  
  coherence_time - MATRIX(nst, ntraj) - for each state is how long has that state resided in a coherence evolution
  coherence_interval - MATRIX(1, ntraj) - for each trajectory - what is the longest coherence time in the system
  decoherence_event_option : 0 - compare the coherence time counter with the decoherence time (simplified DISH)
                             1 - compare the coherence time counter with the time drawn from the exponential distribution
                                 with the parameter lambda = 1/decoherence time - this distribution corresponds to 
                                 the statistics of wait times between the Poisson-distributed events (decoherence)
                                 This is what the original DISH meant to do 
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
      double tau = 0.0;
      if(decoherence_event_option==0){  tau = coherence_interval.get(i, traj); }
      else if(decoherence_event_option==1){  tau = rnd.exponential(1.0/coherence_interval.get(i, traj)); }
        
      if(coherence_time.get(i, traj) >= tau ) {     which_decohere.push_back(i);    }
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

vector<int> decoherence_event(MATRIX& coherence_time, MATRIX& coherence_interval, Random& rnd){

  return decoherence_event(coherence_time, coherence_interval, 0, rnd);

}


vector<int> dish(dyn_control_params& prms,
       MATRIX& q, MATRIX& p,  MATRIX& invM, CMATRIX& Coeff, /*vector<CMATRIX>& projectors, */
       nHamiltonian& ham, vector<int>& act_states, MATRIX& coherence_time, 
       vector<MATRIX>& decoherence_rates, Random& rnd){

    int collapse_option = 0;

    int i,j, traj;
    int nst = Coeff.n_rows; 
    int ntraj = Coeff.n_cols;

    vector<int> old_states(ntraj,-1);
    vector<int> new_states(ntraj,-1);
    vector<int> proposed_states(ntraj,-1); // working variable
    vector<int> which_trajectories(1,-1);

    /// By default, the proposed states are assumed to be the current ones
    vector<int> final_states(act_states);  // this will be the result


    /// Update coherence intervals 
    MATRIX coherence_interval(nst, ntraj); 
    coherence_interval = coherence_intervals(Coeff, decoherence_rates);
 
    /// Determine which states may experience decoherence event
    /// If the decohered_states[traj] == -1, this means no basis states on the trajectory traj
    /// have experienced the decoherence event, othervise the variable will contain an index of 
    /// such state for each trajectory
    vector<int> decohered_states( decoherence_event(coherence_time, coherence_interval, prms.dish_decoherence_event_option, rnd) );

    //cout<<"======= In dish... ==========\n";

    for(traj=0; traj < ntraj; traj++){

      int istate = decohered_states[traj];

      //cout<<"   == Trajectory "<<traj<<"  decohered state: "<<istate<<" ===\n";

      /// Exclude the situation when no decoherence event occurs (-1)
      /// in those cases we just continue the coherent evolution
      if(istate>-1){

        /// No matter what happens with the wavefunctions or whether the hops are accepted,
        /// the cohrence interval is reset for the decohered state, since the state 
        /// has experienced a decoherence event
        coherence_time.set(istate, traj, 0.0);


        /// For both situations below, we'll be making some decisions stochastically, so the following
        /// variables will be needed:
        double prob = (std::conj(Coeff.get(istate, traj)) * Coeff.get(istate, traj) ).real();
        double ksi = rnd.uniform(0.0, 1.0);


        /// Now handle 2 situations: if the decohered state is the active one or if it is not
        ///
        /// Situation 1: Decohered state is the active state
        if(istate == act_states[traj]){

            //cout<<"=== Place 1: decohered state is the active one\n";

            if(ksi<=prob){   
              //cout<<" ... collapsing onto state "<<istate<<endl;
              /// Collapse the wavefunction onto the current active state, stay on that state - no hops
              collapse(Coeff, traj, act_states[traj], collapse_option); 
            }
            else{

              /// Going to try to project out the current state and try hopping onto 
              /// some other state

              vector<int> possible_outcomes;
              possible_outcomes = where_can_we_hop(traj, prms, q, p, invM, Coeff, /*projectors,*/ ham, act_states, rnd);
             
              /// Transitions to some of the states are possible: then, we'll just need to determine where
              int num_allowed = possible_outcomes.size();

              //cout<<" ... can collapse onto: ";
              //for(j=0;j<num_allowed;j++){  cout<<possible_outcomes[j]<<"  "; }
              //cout<<endl;

              if(num_allowed > 0){

                /// Hop into one of the allowed states:
                MATRIX hop_probabilities(1, num_allowed);
 
                double norm = 0.0;
                for(int j=0;j<num_allowed;j++){
                    int jstate = possible_outcomes[j];
                    double pj = (std::conj(Coeff.get(jstate, traj)) * Coeff.get(jstate, traj) ).real();
                    norm += pj;
                    hop_probabilities.set(0, j, pj );
                }
                hop_probabilities *= 1.0/norm;


                double ksi2 = rnd.uniform(0.0, 1.0);
                int selected_state_index = hop(0, hop_probabilities, ksi2);

                cout<<"  selected_state_index = "<<selected_state_index<<"  value = "<<possible_outcomes[selected_state_index]<<endl;

                /// We hop to this state
                final_states[traj] = possible_outcomes[selected_state_index];

                /// We project out the current active state
                project_out(Coeff, traj, istate); 

              }/// num_allowed > 0
              /// There are no states to where the transitions are possible => do nothing
              else{ ;; } ///  Do nothing

            }

        }

        /// Situation 2: Decohered state is an inactive state
        else{
            //cout<<"=== Place 2: decohered state is one of the inactive states\n";

            if(ksi<=prob){   
              /// The decohered state becomes the proposed state
              proposed_states[traj] = istate;
              which_trajectories[0] = traj;

              /// Decide if we can accept the transitions, the function below only checks the hopping for a single trajectory `traj`
              /// other elements of the input and output vector<int> variables (old_states, new_states, proposed_states) are irrelevant
              /// the variable `which_trajectories` instructs to handle only the current trajectory
              old_states = act_states;
              new_states = accept_hops(prms, q, p, invM, Coeff, /*projectors,*/ ham, proposed_states, old_states, rnd, which_trajectories);              

           
              //cout<<"proposed_state = "<<proposed_states[traj]<<"  accepted hop "<<new_states[traj]<<endl;

              /// The proposed transition can be accepted => collapse onto this state, hop onto this state
              if(new_states[traj] == proposed_states[traj]){

                /// Successfull hop - collapse onto this new state
                collapse(Coeff, traj, new_states[traj], collapse_option); 

                /// Mark the hop
                final_states[traj] = new_states[traj];

              }
              else{  ;; }  /// Do nothing

            }
            else{
              /// Project out the state from the superposition, no hops
              project_out(Coeff, traj, istate); 
            }

        }/// decohered state is the inactive state


      }// istate > -1
    
    }// for traj

    return final_states;

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
          project_out(Coeff, traj, decohered_states[traj]); 
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

