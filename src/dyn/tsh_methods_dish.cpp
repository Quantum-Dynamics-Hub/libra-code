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
  \file tsh_methods_dish.cpp
  \brief The file implements Decoherence Induced Surface Hopping (DISH) algorithm
    
*/

#include "Surface_Hopping.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


MATRIX coherence_intervals(CMATRIX& Coeff, MATRIX& rates){
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

