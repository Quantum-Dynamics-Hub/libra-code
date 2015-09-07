/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "Dynamics_Nuclear.h"
#include "Energy_and_Forces.h"


namespace libdyn{



void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,int opt){

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


}// propagate_nuclear





}// namespace libdyn
