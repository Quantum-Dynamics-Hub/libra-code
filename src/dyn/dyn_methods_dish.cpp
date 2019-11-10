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

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{


void dish(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
          nHamiltonian& ham, bp::object py_funct, bp::dict params, bp::dict dyn_params, Random& rnd){
/**
  \brief Decoherence-induced surface hopping (DISH)
  Reference: Jaeger, H. M.; Fischer, S.; Prezhdo, O. V. Decoherence-Induced Surface Hopping. J. Chem. Phys. 2012, 137, 22A545.

  This function actually combines 

  \param[in,out] q - coordinates of nuclei
  \param[in,out] Cadi - amplitudes of the adiabatic states
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

  dyn_control_params prms;
  prms.set_parameters(dyn_params);

/*


  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int nst = Cadi.n_rows;    
  int traj, dof;


  const double kb = 3.166811429e-6; // Hartree/K   
  int i,j; 
  int has_decoherence = 0; /// set to 1 if we have already encountered a decoherence event

  for(i=0;i<nst && !has_decoherence;i++){

    /// The state i has evolved coherently for longer than the coherence interval
    /// so it has to experience a decoherence event 
    if(t_m.get(i) >= tau_m.get(i) ) { 


    for(traj=0; traj<ntraj; traj++){
      prev_ham_dia[traj] = ham.children[traj]->get_ham_dia().real();  
    }


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
*/

} // dish




}// namespace libdyn
}// liblibra

