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
  \file dyn_decoherence_methods.cpp
  \brief The file implements various decoherence correction methods
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


CMATRIX sdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates){
    /**
    \brief Generic framework of Simplified Decay of Mixing (SDM) method of 
    Granucci, G.; Persico, M. J. Chem. Phys. 2007, 126, 134114
    
    \param[in]       Coeff [ CMATRIX(nadi, 1) ] An object containig electronic DOFs. 
    \param[in]          dt [ float ] The integration timestep. Units = a.u. of time
    \param[in]      act_st [ integer ] The active state index
    \param[in]      decoh_rates [ MATRIX ] The matrix of decoherence (pure dephasing) rates between all pairs of states

    The function returns:
    C [ CMATRIX ] - the updated state of the electronic DOF, in the same data type as the input

    */

    double sclf;


    CMATRIX C(Coeff);


    // Population of the active state
    double p_aa_old = (std::conj(C.get(act_st)) * C.get(act_st)).real();

    if(p_aa_old>1.0){
      // Comment this place if you want to allow inprecise integration
      // (where the total norm may exceeed 1.0), as is the case for too large dt
      // for some algorithms
      cout<<"=== Place 1 =====\n";
      cout<<"Error in CMATRIX msdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
      cout<<"active state is larger than 1: p_aa_old = "<< p_aa_old << endl;
      cout<<"C = \n"; C.show_matrix();
      cout<<"act_st = "<<act_st<<endl;
      cout<<"Coeff = \n"; Coeff.show_matrix();
      cout<<"decoh_rates = \n"; decoh_rates.show_matrix();
      cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
      exit(0);
    }


    if(p_aa_old>0.0){

      // First - update all the coefficients for the non-active states        
      int N = Coeff.n_elts;

      double inact_st_pop = 0.0; // population of the inactive states after rescaling

      for(int i=0; i<N; i++){
        if(i != act_st){
          double itau = decoh_rates.get(i, act_st); 
          sclf = exp(-dt*itau);
          C.scale(i, 0, sclf);

          inact_st_pop += (std::conj(C.get(i)) * C.get(i)).real();
        }
      }

      if(inact_st_pop>1.0){
        cout<<"=== Place 2 =====\n";
        cout<<"Error in CMATRIX sdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
        cout<<"Total population of inactive states after rescaling is larger than 1: inact_st_pop = "<<inact_st_pop<<endl;
        cout<<"C = \n"; C.show_matrix();
        cout<<"act_st = "<<act_st<<endl;
        cout<<"Coeff = \n"; Coeff.show_matrix();
        cout<<"decoh_rates = \n"; decoh_rates.show_matrix();
        cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
        exit(0);
      }

      double p_aa_new = 1.0 - inact_st_pop;

     
      if(p_aa_new<0.0){
        cout<<"=== Place 3 =====\n";
        cout<<"Error in CMATRIX sdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
        cout<<"new population of the active state "<< p_aa_new <<" is negative...\n";
        cout<<"inact_st_pop = "<<inact_st_pop<<endl;
        cout<<"p_aa_old = "<<p_aa_old<<endl;
        cout<<"C = \n"; C.show_matrix();
        cout<<"act_st = "<<act_st<<endl;
        cout<<"Coeff = \n"; Coeff.show_matrix();
        cout<<"decoh_rates = \n"; decoh_rates.show_matrix();     
        cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
        exit(0);
      }

      sclf = sqrt( p_aa_new / p_aa_old );  // scaling factor for the active state
        
      // Rescale the active state
      C.scale(act_st, 0, sclf);

    }// if p_aa_old > 0.0

    double new_norm = (C.H() * C).get(0,0).real();
    //cout<<"new_norm = "<<new_norm<<endl;


    if(fabs(new_norm-1.0)>0.1){
      cout<<"=== Place 4 =====\n";
      cout<<"Error in CMATRIX sdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates):\n";
    //  cout<<"new population of the active state "<< p_aa_new <<" is negative...\n";
    //  cout<<"inact_st_pop = "<<inact_st_pop<<endl;
      cout<<"p_aa_old = "<<p_aa_old<<endl;
      cout<<"C = \n"; C.show_matrix();
      cout<<"act_st = "<<act_st<<endl;
      cout<<"Coeff = \n"; Coeff.show_matrix();
      cout<<"decoh_rates = \n"; decoh_rates.show_matrix();     
      cout<<"initial total pop = "<<(Coeff.H() * Coeff).get(0,0).real();
      exit(0);
    }


    return C;

}


CMATRIX sdm(CMATRIX& Coeff, double dt, vector<int>& act_st, vector<MATRIX>& decoh_rates){
    /**
    \brief The generic framework of the Simplified Decay of Mixing (SDM) method of
    Granucci, G.; Persico, M. J. Chem. Phys. 2007, 126, 134114)

    This is a version for multiple trajectories
    
    \param[in]       Coeff [ CMATRIX(nadi, ntraj) ] An object containig electronic DOFs. 
    \param[in]          dt [ float ] The integration timestep. Units = a.u. of time
    \param[in]      act_st [ integer ] The active state index
    \param[in]      decoh_rates [ MATRIX ] The matrix of decoherence (pure dephasing) rates between all pairs of states

    The function returns:
    # C [ CMATRIX ] - the updated state of the electronic DOF, in the same data type as the input

    */

  int nadi = Coeff.n_rows;
  int ntraj = Coeff.n_cols;
  int i, traj;


  vector<int> stenc_x(nadi, 0); for(i=0;i<nadi;i++){  stenc_x[i] = i; }
  vector<int> stenc_y(1, 0); 

  CMATRIX coeff(nadi, 1);
  CMATRIX res(nadi, ntraj);
 
  for(traj=0; traj<ntraj; traj++){

    stenc_y[0] = traj;
    pop_submatrix(Coeff, coeff, stenc_x, stenc_y);
    coeff = sdm(coeff, dt, act_st[traj], decoh_rates[traj]);
    push_submatrix(res, coeff, stenc_x, stenc_y);

  }// for traj
  return res;

}



void project_out(CMATRIX& Coeff, int traj, int i){
    /** Projects the state `i` out of a coherent superposition of states

    Args:
        Coeff ( CMATRIX(nstates, ntraj) ): amplitudes of the electronic states in
            a coherent superposition from which we will project a state out
        traj ( int ): The index of the trajectory (column) on which we will operate
        i ( int ): The index of the state to be projected out

    Returns:
        None: but changes the input variable `Coeff`

   */ 

    int nstates = Coeff.n_rows;
    
    complex<double> ci; ci = Coeff.get(i,traj);
    double pi; pi = (std::conj(ci) * ci).real();
    double nrm = 1.0 - pi;

    if(nrm<=0.0){ nrm = 0.0; }
    if(nrm>0.0){  nrm = 1.0/sqrt(nrm); }

    Coeff.scale(-1, traj, nrm);
    Coeff.set(i, traj, complex<double>(0.0, 0.0));
}


void collapse(CMATRIX& Coeff, int traj, int i, int collapse_option){
    /** Collapse the wfc but such that to preserve the phase!

    Args:
        Coeff ( CMATRIX(nstates, ntraj) ): amplitudes of the electronic states in
            a coherent superposition which we are going to collapse
        traj ( int ): The index of the trajectory (column) on which we will operate
        i ( int ): The index of the state onto which the supeposition will be collapsed
        collapse_option ( int ): the way to collapse onto given state:
        
          - 0: by rescaling the magnitude of the amplitude vector elements, but preserving "phase"
          - 1: by resetting the amplitudes to 1.0+0.0j. This option changes phase

    Returns:
        None: but changes the input variable `Coeff`

    */ 

    complex<double> ci; ci = Coeff.get(i,traj);
    double pi; pi = (std::conj(ci) * ci).real();

    Coeff.scale(-1, traj, 0.0);

    if(collapse_option==0){
      if(pi>0.0){   Coeff.set(i, traj, ci/sqrt(pi)); }     
      else{  Coeff.set(i, traj, complex<double>(1.0, 0.0) ); }
    }
    else if(collapse_option==1){
      Coeff.set(i, traj, complex<double>(1.0, 0.0) );
    }

}



void instantaneous_decoherence(CMATRIX& Coeff, 
   vector<int>& accepted_states, vector<int>& proposed_states, vector<int>& initial_states,
   int instantaneous_decoherence_variant, int collapse_option){

/**
  This function implements the instantaneous decoherence approach of Nelson et al.:

  Nelson, T.; Fernandez-Alberti, S.; Roitberg, A. E.; Tretiak, S. J. Chem. Phys. 2013, 138, 224111

  Two options of the algorithm are available:

  ID-S: wavefunction amplitudes are collapsed only during the successful hops (onto the new state)

  ID-A: wavefunction amplitudes are collapsed at every attempted hop
      to the new state, if successful
      to the old state, if now

  There collapsing options are controlled by the parameter instantaneous_decoherence_variant

   0 - ID-S
   1 - ID-A
   2 - ID-C - consistent ID - an experimental algo

   In the "consistent" version, we lift the condition that the accepted/proposed states must be different from
   the starting state - this option addresses a philosophical question - what if we say that no hops 
   is equivalent to the hop onto the current state? why should the "hopping" into the original state 
   by treated differently from the "actual hopping" to another state? So, in this version we 
   collapse onto the accpeted states, no matter if they are the result of a successfull or frustrated hop.

  The mechanism of the collapse event itself is controlled by the collapse_option parameter
*/
  int traj;
  int ntraj = Coeff.n_cols;

  if(accepted_states.size()!=ntraj){
    cout<<"ERROR in ids: the sizes of the input variables Coeff and accepted_states are inconsistent\n";
    cout<<"Coeff.num_of_cols = = "<<ntraj<<"\n";
    cout<<"accepted_states.size() = "<<accepted_states.size()<<"\n";
    cout<<"exiting...\n";
    exit(0);
  }
  if(proposed_states.size()!=ntraj){
    cout<<"ERROR in ids: the sizes of the input variables Coeff and proposed_states are inconsistent\n";
    cout<<"Coeff.num_of_cols = = "<<ntraj<<"\n";
    cout<<"proposed_states.size() = "<<proposed_states.size()<<"\n";
    cout<<"exiting...\n";
    exit(0);
  }
  if(initial_states.size()!=ntraj){
    cout<<"ERROR in ids: the sizes of the input variables Coeff and initial_states are inconsistent\n";
    cout<<"Coeff.num_of_cols = = "<<ntraj<<"\n";
    cout<<"initial_states.size() = "<<initial_states.size()<<"\n";
    cout<<"exiting...\n";
    exit(0);
  }


  // ID-S
  if(instantaneous_decoherence_variant==0){

    for(traj = 0; traj < ntraj; traj++){

      if(accepted_states[traj] != initial_states[traj]){
        collapse(Coeff, traj, accepted_states[traj], collapse_option);
      }
    }// traj

  }// ID-S

  // ID-A
  else if(instantaneous_decoherence_variant==1){

    for(traj = 0; traj < ntraj; traj++){

      if(proposed_states[traj] != initial_states[traj]){
      // Only apply ID-A, if the proposed states are different from the original ones

        if(accepted_states[traj] == proposed_states[traj]){
        // Proposed hop is successful - collapse onto newly accepted state
          collapse(Coeff, traj, accepted_states[traj], collapse_option);
        }
        else{
        // Proposed hop is not successful - collapse onto the original state
          collapse(Coeff, traj, initial_states[traj], collapse_option);
        }
      }
    }// traj

  }// ID-A


  // ID-C
  else if(instantaneous_decoherence_variant==2){

    for(traj = 0; traj < ntraj; traj++){

      collapse(Coeff, traj, accepted_states[traj], collapse_option);

    }// traj

  }// ID-C



}


}// namespace libdyn
}// liblibra

