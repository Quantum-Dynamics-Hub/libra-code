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
  \file Dynamics_Nuclear.cpp
  \brief The file implements the functions for nuclear (classical) dynamics
    
*/

#include "Dynamics_Nuclear.h"
#include "Energy_and_Forces.h"


/// liblibra namespace
namespace liblibra{

/// libdyn namespace 
namespace libdyn{



void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,int opt){
/**
  \brief One step of Velocity verlet algorithm for nuclear DOF
  \param[in] Integration time step
  \param[in,out] mol Describes the nuclear DOF. Changes during the integration.
  \param[in] el Describes electronic DOF. Does not change during the integration.
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

*/


  int i;
  vector<double> v_(mol->nnucl); 

  // Propagate momenta and positions of all nuclei
  // exp(iL_q*dt)*exp(iL_p*dt/2)
  mol->propagate_p(0.5*dt);
  mol->propagate_q(dt);

  // Update Hamiltonian
  ham->set_q(mol->q);
  for(i=0;i<mol->nnucl;i++){ v_[i] = mol->p[i]/mol->mass[i]; }  ham->set_v(v_);
  ham->compute();  

  // Update forces and potential energy
  compute_potential_energy(mol, el, ham, opt);
  compute_forces(mol, el, ham, opt);


  // Propagate momenta of all nucleii
  // operator exp(iL_p*dt/2)
  mol->propagate_p(0.5*dt);


//  ham->set_q(mol->q);
  for(i=0;i<mol->nnucl;i++){ v_[i] = mol->p[i]/mol->mass[i]; }  ham->set_v(v_);
//  ham->compute();  



}// propagate_nuclear


void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,Thermostat* therm, int opt){
/**
  \brief One step of Velocity verlet algorithm coupled to nuclear Thermostat - for NVT calculations
  !!!  IMPORTANT: Presently, this is only a template - the thermostat part is only outlined - need to test this first !!!

  \param[in] Integration time step
  \param[in,out] mol Describes the nuclear DOF. Changes during the integration.
  \param[in] el Describes electronic DOF. Does not change during the integration.
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) -its internal variables
  are changed during the compuations
  \param[in] therm The pointer to a Thermostat object, which contains the information about thermal bath state
  \param[in] opt Option for selecting the way to describe the electron-nuclear interaction: = 0 - Ehrenfest (mean-field, MF), =1 -
  (fewest switched surface hopping, FSSH)

*/

  int i;

  // Propagate momenta and positions of all nuclei
//  mol->P[i] *= exp(-therm->gamma*0.25*dt);
//  mol->P[i] += therm->ksi() * 0.5*dt;
//  mol->P[i] *= exp(-therm->gamma*0.25*dt);


  // exp(iL_q*dt)*exp(iL_p*dt/2)
  mol->propagate_p(0.5*dt);
  mol->propagate_q(dt);

  // Update Hamiltonian
  ham->set_q(mol->q);  ham->compute();  

  // Update forces and potential energy
  compute_potential_energy(mol, el, ham, opt);
  compute_forces(mol, el, ham, opt);


  // Propagate momenta of all nucleii
  // operator exp(iL_p*dt/2)
  mol->propagate_p(0.5*dt);

  // Propagate momenta and positions of all nuclei
//  mol->P[i] *= exp(-therm->gamma*0.25*dt);
//  mol->P[i] += therm->ksi() * 0.5*dt;
//  mol->P[i] *= exp(-therm->gamma*0.25*dt);



}// propagate_nuclear


}// namespace libdyn
}// liblibra


