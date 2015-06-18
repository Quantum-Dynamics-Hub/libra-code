#include "Dynamics_Nuclear.h"
#include "Energy_and_Forces.h"


namespace libdyn{



void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,int opt){

  int i;

  // Propagate momenta and positions of all nuclei
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
