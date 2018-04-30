/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Surface_Hopping.cpp
  \brief The file implements the functions used in surface hopping methods
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


void tsh(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, MATRIX& states,
         nHamiltonian& ham, bp::object py_funct, bp::object params, int rep){
 
  //============== Electronic propagation ===================
  if(rep==0){  
    ham.compute_nac_dia(p, invM);
    ham.compute_hvib_dia();
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM); 
    ham.compute_hvib_adi();
  }

  propagate_electronic(0.5*dt, C, ham, rep);   

  //============== Nuclear propagation ===================
    
       if(rep==0){  p = p + ham.Ehrenfest_forces_dia(C).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(C).real() * 0.5*dt;  }


  q = q + invM*p*dt;
  ham.compute_diabatic(py_funct, bp::object(q), params);
  ham.compute_adiabatic(1);


       if(rep==0){  p = p + ham.Ehrenfest_forces_dia(C).real() * 0.5*dt;  }
  else if(rep==1){  p = p + ham.Ehrenfest_forces_adi(C).real() * 0.5*dt;  }

  //============== Electronic propagation ===================
  if(rep==0){  
    ham.compute_nac_dia(p, invM);
    ham.compute_hvib_dia();
  }
  else if(rep==1){  
    ham.compute_nac_adi(p, invM); 
    ham.compute_hvib_adi();
  }

  propagate_electronic(0.5*dt, C, ham, rep);   


}






int hop(int initstate, MATRIX& g, double ksi){
/** 
  \brief Attempts a stochastic hop from the initial state "initstate"
  \param[in] initstate The index of the state from which we try to hop out 
  \param[in] g The hopping probabilities matrix (type MATRIX)
  \param[in] ksi A random number that determines the outcome of the "hop" procedure

  Returned value: the index of the state to which we have hopped
*/

  int nstates = g.n_cols;
  double left, right; left = right = 0.0;
  int finstate = initstate;

  for(int i=0;i<nstates;i++){
    if(i==0){left = 0.0; right = g.get(initstate,i); }
    else{  left = right; right = right + g.get(initstate,i); }
 
    if((left<ksi) && (ksi<=right)){  finstate = i;  }
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



int ida(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi){
/**
  \brief Instantaneous decoherence at attempted hops (ID-A)

  \param[in,out] Coeff The matrix of coefficients (amplitudes of the basis excited states). This matrix is modified at every
  decoherence event
  \param[in] old_st The old state before an attempted hop
  \param[in] new_st The new state after an attempted hop (no matter what the algorithm - FSSH, GFSH, MSSH, or something else)
  \param[in] E_old The energy corresponding to the old state
  \param[in] E_new The energy corresponding to the new state
  \param[in] T is the temperature of the system
  \param[in] ksi a random number from a uniform distribution of the [0,1] interval - needed to decide the hop acceptance/decoherence

  The function returns the index of the new state after hop rejection/decoherence criteria
  The function also modifies the amplitudes of the coherent superposition

*/

  const double kb = 3.166811429e-6;  // Boltzmann constant: Hartree/K
  int istate = old_st;               // by default the resulting state is the old state
  
  if(new_st != old_st){  // attempted hop

    // Now apply energy considerations
    double dE = (E_new - E_old);
    double boltz_f = 1.0;

    if(dE>0.0){
      double argg = dE/(kb*T);
      if(argg > 50.0){ boltz_f = 0.0; }
      else{            boltz_f = exp(-argg); }

      if(ksi<boltz_f){
        istate = new_st;  // accepted hop

        // Collapse the wavefunction onto the new state
        Coeff *= 0.0;
        Coeff.set(new_st, 1.0, 0.0);
      }
      else{
        // Unsuccessful hop - collapse wfc back to the original state
        Coeff *= 0.0;
        Coeff.set(old_st, 1.0, 0.0);
      }
    }// dE > 0
  }// new_st != old_st

  return istate;

}


MATRIX coherence_intervals(CMATRIX& Coeff, MATRIX& rates ){
/**
  This function computes the time-dependent (and population-dependent) coherence intervals
  (the time after which different states should experience a decoherence event)
  as described by Eq. 11 in:
  Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  1/tau_i  (t) =  sum_(j!=i)^nstates {  rho_ii(t) * rate_ij }


  \param[in] Coeff Amplitudes of the electronic states
  \param[in] rates A matrix containing the decoherence rates (inverse of the
  decoherence time for each given pair of states)

  Returns: A matrix of the coherence intervals for each state

*/
  int nstates = Coeff.n_rows; 

  CMATRIX* denmat; denmat = new CMATRIX(nstates, nstates);   
  *denmat = (Coeff * Coeff.H() ).conj();

  MATRIX tau_m(nstates, 1); tau_m *= 0.0;
  

  for(int i=0;i<nstates;i++){

    double summ = 0.0;
    for(int j=0;j<nstates;j++){

      if(j!=i){
        summ += denmat->get(j,j).real() * rates.get(i,j); 
      }// if

    }// for j

    if(summ>0.0){   tau_m.set(i, 1.0/summ); }
    else        {   tau_m.set(i, 1.0e+100); } // infinite coherence interval
    
     
  }// for i

  delete denmat;

  return tau_m;
}



int dish(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib, int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2){
/**
  \brief Decoherence-induced surface hopping (DISH)

  Reference: Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  \param[in,out] el Electronic object, containing the info about electronic amplitudes
  \param[in,out] t_m A matrix N x 1 of the times each state resides in a coherence interval
   (since the last decoherence event)
  \param[in] tau_m A matrix N x 1 of the coherence intervals for each electronic state
  \param[in] Hvib The matrix of vibronic Hamiltonian - we need the energies
  \param[in] use_boltz_flag  if set to 1, the hopping probabilities will be re-scaled by the 
  Boltzmann factor. This is need in neglect of back-reaction approximation (NBRA) calculations, 
  when no 
  \param[in] Ekin Kinetic energy of nuclei
  \param[in] T is the temperature of the system
  \param[in] ksi1 a is random number from a uniform distribution on the [0,1] interval
  \param[in] ksi2 a is random number from a uniform distribution on the [0,1] interval

  The function modifies  el and t_m variables
  Returns: the index of the electronic state after potential hop
*/

  const double kb = 3.166811429e-6; // Hartree/K   
  int i,j; 
  int has_decoherence = 0; /// set to 1 if we have already encountered a decoherence event

  for(i=0;i<el.nstates && !has_decoherence;i++){

    /// The state i has evolved coherently for longer than the coherence interval
    /// so it has to experience a decoherence event 
    if(t_m.get(i) >= tau_m.get(i) ) { 

      /// There are essentially two outcomes when the decoherence takes place:
      /// One: we collapse the wavefunction onto the state i with the probability 
      /// given by the population of that state

      if(ksi1 < el.rho(i,i).real()){ 

        /// Now, lets determine if the hop is possible based on the energy conservation
        /// considerations 



        int can_hop = 0;

        double E_i = Hvib.get(el.istate, el.istate).real();/// initial potential energy
        double E_f = Hvib.get(i,i).real();                 /// proposed potential energy
        double dE = E_f - E_i;

        /// In leu of hop rejection use Boltzmann factors: for NBRA simulations
        if(use_boltz_flag==1){
          double bf = 1.0;
          if(dE>0){  bf = exp(-dE/(kb*T)); }  /// hop to higher energy state is difficult
          if(ksi2<bf){ can_hop = 1; }  /// hop is allowed thermally

        }
        /// Regular energy-conservation based criterion
        else{   
          /// Predicted final kinetic energy is positive - hop is possible
          if(Ekin - dE > 0.0){  can_hop = 1; }
        }
          
        /// Now, decide about decoherence         
        if(can_hop){     el.collapse(i, 1);  } /// Here is the actuall collapse
        else{ el.project_out(i); }

      }

      /// Second: project the system out of that state otherwise
      else{  el.project_out(i);   }


      /// Reset the time axis for state i (only for this state)
      /// other states still reside in a coherent superposition
      t_m.set(i, 0.0);

      /// Set the flag that we have attempted a decoherence event
      /// so we done with DISH at this point in time
      has_decoherence = 1;


    }// t_m[i]>=1.0/tau_m
      
  }// for i

  return el.istate;

} // dish


int dish(Electronic& el, Nuclear& mol, Hamiltonian& ham, MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2){
/**
  \brief Decoherence-induced surface hopping (DISH) - overloaded version

  Reference: Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  \param[in,out] el Electronic object, containing the info about electronic amplitudes
  \param[in,out] t_m A matrix N x 1 of the times each state resides in a coherence interval
   (since the last decoherence event)
  \param[in] tau_m A matrix N x 1 of the coherence intervals for each electronic state
  \param[in] Hvib The matrix of vibronic Hamiltonian - we need the energies
  \param[in] use_boltz_flag  if set to 1, the hopping probabilities will be re-scaled by the 
  Boltzmann factor. This is need in neglect of back-reaction approximation (NBRA) calculations, 
  when no 
  \param[in] Ekin Kinetic energy of nuclei
  \param[in] T is the temperature of the system
  \param[in] ksi1 a is random number from a uniform distribution on the [0,1] interval
  \param[in] ksi2 a is random number from a uniform distribution on the [0,1] interval

  The function modifies  el and t_m variables
  Returns: the index of the electronic state after potential hop
*/

  const double kb = 3.166811429e-6; // Hartree/K   
  int i,j; 
  int has_decoherence = 0; /// set to 1 if we have already encountered a decoherence event

  for(i=0;i<el.nstates && !has_decoherence;i++){

    /// The state i has evolved coherently for longer than the coherence interval
    /// so it has to experience a decoherence event 
    if(t_m.get(i) >= tau_m.get(i) ) { 

      /// There are essentially two outcomes when the decoherence takes place:
      /// One: we collapse the wavefunction onto the state i with the probability 
      /// given by the population of that state

      if(ksi1 < el.rho(i,i).real()){ 

        /// Now, lets determine if the hop is possible based on the energy conservation
        /// considerations 



        int can_hop = 0;

        double E_i = ham.Hvib(el.istate, el.istate).real();/// initial potential energy
        double E_f = ham.Hvib(i,i).real();                 /// proposed potential energy
        double dE = E_f - E_i;

        /// In leu of hop rejection use Boltzmann factors: for NBRA simulations
        if(use_boltz_flag==1){
          double bf = 1.0;
          if(dE>0){  bf = exp(-dE/(kb*T)); }  /// hop to higher energy state is difficult
          if(ksi2<bf){ can_hop = 1; }  /// hop is allowed thermally

        }
        /// Regular energy-conservation based criterion
        else{   
          /// Predicted final kinetic energy is positive - hop is possible
          double Ekin = compute_kinetic_energy(mol); // initial kinetic energy
          if(Ekin - dE > 0.0){  can_hop = 1; }
        }
          
        /// Now, decide about decoherence         
        if(can_hop){     el.collapse(i, 1);  } /// Here is the actuall collapse
        else{ el.project_out(i); }

      }

      /// Second: project the system out of that state otherwise
      else{  el.project_out(i);   }


      /// Reset the time axis for state i (only for this state)
      /// other states still reside in a coherent superposition
      t_m.set(i, 0.0);

      /// Set the flag that we have attempted a decoherence event
      /// so we done with DISH at this point in time
      has_decoherence = 1;


    }// t_m[i]>=1.0/tau_m
      
  }// for i

  return el.istate;

} // dish




}// namespace libdyn
}// liblibra

