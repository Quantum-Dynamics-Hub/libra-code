/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Verlet_NHC.cpp
  \brief The file implements the velocity Verlet algorithm for nuclear dynamics integration coupled
  to a Nose-Hoover Chain thermostat    
*/

#include "Dynamics.h"
#include "Energy_and_Forces.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;


/// libdyn namespace 
namespace libdyn{

using namespace libthermostat;


void Verlet0_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, Thermostat& therm){
/**
  \brief One step of velocity Verlet algorithm for nuclear DOF

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in,out] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) - its internal variables
  (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] lvl The level of Hamiltonians to do the evaluation at

*/

  int ndof = q.n_rows;
  int dof;

  int ham_rep = 0; // default -- assume the Hamiltonian is first computed in the diabatic representation
                   // and then will be transformed to the adiabatic in this function. 

  std::string key;
  boost::python::dict d = (boost::python::dict)params;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);
    if(key=="ham_rep") { ham_rep = extract<int>(d.values()[i]);   }
  }



  CMATRIX Cadi(1,1); Cadi.set(0,0,1.0,0.0);

  
  p *= therm.vel_scale(0.5*dt);     
  p = p + ham.forces_adi(Cadi).real() * 0.5*dt;

  // For efficiency, switch to the element-wise multiplication
  for(dof=0; dof<ndof; dof++){  
    q.add(dof, 0,  invM.get(dof,0) * p.get(dof,0) * dt ); 
  }

  if(ham_rep==0){
    ham.compute_diabatic(py_funct, bp::object(q), params, 0);
    ham.compute_adiabatic(1, 0);
  }
  else if(ham_rep==1){
    ham.compute_adiabatic(py_funct, bp::object(q), params, 0);
  }

  double ekin = compute_kinetic_energy(p, invM);
  therm.propagate_nhc(dt, ekin, 0.0, 0.0);

  p = p + ham.forces_adi(Cadi).real() * 0.5*dt;
  p *= therm.vel_scale(0.5*dt);     

}// Verlet0



void Verlet1_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, 
                 int entanglement_opt, vector<Thermostat>& therm){
/**
  \brief One step of velocity Verlet algorithm for nuclear DOF

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in,out] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) - its internal variables
  (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] lvl The level of Hamiltonians to do the evaluation at
  \param[in] entanglement_opt - a selector of a method to couple the trajectories in this ensemble:
             0 - no coupling, 1 - ETHD, 2 - RPMD

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int traj, dof;

  int ham_rep = 0; // default -- assume the Hamiltonian is first computed in the diabatic representation
                   // and then will be transformed to the adiabatic in this function. 


  //============= Extract optional parameters: needed for some execution scenarios =============
  double ETHD3_alpha = 1.0;
  std::string key;
  boost::python::dict d = (boost::python::dict)params;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);
    if(key=="ETHD3_alpha") { ETHD3_alpha = extract<double>(d.values()[i]);   }
    if(key=="ham_rep") { ham_rep = extract<int>(d.values()[i]);   }
  }


  vector<int> t1(ndof, 0); for(int i=0;i<ndof;i++){  t1[i] = i; }
  vector<int> t2(1,0);
  vector<int> t3(2,0);

  CMATRIX Cadi(1,1); Cadi.set(0,0,1.0,0.0);
  MATRIX F(ndof, ntraj);
  MATRIX f(ndof, 1);


  for(traj=0; traj<ntraj; traj++){
    p.scale(-1, traj, therm[traj].vel_scale(0.5*dt));
  }
  
  for(traj=0; traj<ntraj; traj++){
    t2[0] = traj;  t3[1] = traj;
    f = ham.forces_adi(Cadi, t3).real();
    push_submatrix(F, f, t1, t2);
  }


  p = p + F * 0.5*dt;

  // For efficiency, switch to the element-wise multiplication
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  

      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * dt ); 

    }
  }

  if(ham_rep==0){
    ham.compute_diabatic(py_funct, bp::object(q), params, 1);
    ham.compute_adiabatic(1, 1);
  }
  else if(ham_rep==1){
    ham.compute_adiabatic(py_funct, bp::object(q), params, 1);
  }


  if(entanglement_opt==0){    /* Nothing to do */   }
  else if(entanglement_opt==1){   ham.add_ethd_adi(q, invM, 1);  }
  else if(entanglement_opt==2){   ham.add_ethd3_adi(q, invM, ETHD3_alpha, 1);  }
  else{
    cout<<"ERROR in Verlet1: The entanglement option = "<<entanglement_opt<<" is not avaialable\n";
    exit(0);
  }

  for(traj=0; traj<ntraj; traj++){
    t2[0] = traj; 
    push_submatrix(p, f, t1, t2);
    double ekin = compute_kinetic_energy(f, invM);
    therm[traj].propagate_nhc(dt, ekin, 0.0, 0.0);
  }


  for(traj=0; traj<ntraj; traj++){
    t2[0] = traj;  t3[1] = traj;
    f = ham.forces_adi(Cadi, t3).real();
    push_submatrix(F, f, t1, t2);
  }

  p = p + F * 0.5*dt;

  for(traj=0; traj<ntraj; traj++){
    p.scale(-1, traj, therm[traj].vel_scale(0.5*dt));
  }




}// Verlet1


void Verlet1_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, 
                 bp::object params, vector<Thermostat>& therm){
/**
  \brief Regular velocity Verlet for an ensemble of (uncoupled trajectories)
*/

  Verlet1_nvt(dt, q, p, invM, ham, py_funct, params, 0, therm);


}// Verlet1





void Verlet1_nvt(double dt, MATRIX& q, MATRIX& P, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, 
                 int entanglement_opt, Thermostat& therm){
/**
  \brief One step of velocity Verlet algorithm for nuclear DOF

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in,out] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) - its internal variables
  (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] lvl The level of Hamiltonians to do the evaluation at
  \param[in] entanglement_opt - a selector of a method to couple the trajectories in this ensemble:
             0 - no coupling, 1 - ETHD, 2 - RPMD

*/

  int ndof = q.n_rows;
  int ntraj = q.n_cols;
  int traj, dof;

  //============= Extract optional parameters: needed for some execution scenarios =============
  double ETHD3_alpha = 1.0;
  std::string key;
  boost::python::dict d = (boost::python::dict)params;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);
    if(key=="ETHD3_alpha") { ETHD3_alpha = extract<double>(d.values()[i]);   }
  }


  vector<int> t1(ndof, 0); for(int i=0;i<ndof;i++){  t1[i] = i; }
  vector<int> t2(1,0);
  vector<int> t3(2,0);

  CMATRIX Cadi(1,1); Cadi.set(0,0,1.0,0.0);
  MATRIX F(ndof, ntraj);
  MATRIX f(ndof, 1);
  MATRIX p(ndof, 1);


  P = P * therm.vel_scale(0.5*dt);
  
  for(traj=0; traj<ntraj; traj++){
    t2[0] = traj;  t3[1] = traj;
    f = ham.forces_adi(Cadi, t3).real();
    push_submatrix(F, f, t1, t2);
  }


  P = P + F * 0.5*dt;

  // For efficiency, switch to the element-wise multiplication
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  

      q.add(dof, traj,  invM.get(dof,0) * P.get(dof,traj) * dt ); 

    }
  }

  ham.compute_diabatic(py_funct, bp::object(q), params, 1);
  ham.compute_adiabatic(1, 1);

  if(entanglement_opt==0){    /* Nothing to do */   }
  else if(entanglement_opt==1){   ham.add_ethd_adi(q, invM, 1);  }
  else if(entanglement_opt==2){   ham.add_ethd3_adi(q, invM, ETHD3_alpha, 1);  }
  else{
    cout<<"ERROR in Verlet1: The entanglement option = "<<entanglement_opt<<" is not avaialable\n";
    exit(0);
  }


  double ekin = 0.0;
  for(traj=0; traj<ntraj; traj++){
    t2[0] = traj; 
    pop_submatrix(P, p, t1, t2);
    ekin += compute_kinetic_energy(p, invM);
  }
  therm.propagate_nhc(dt, ekin, 0.0, 0.0);



  for(traj=0; traj<ntraj; traj++){
    t2[0] = traj;  t3[1] = traj;
    f = ham.forces_adi(Cadi, t3).real();
    push_submatrix(F, f, t1, t2);
  }

  P = P + F * 0.5*dt;
  P = P * therm.vel_scale(0.5*dt);




}// Verlet1


void Verlet1_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, 
                 bp::object params, Thermostat& therm){
/**
  \brief Regular velocity Verlet for an ensemble of (uncoupled trajectories)
*/

  Verlet1_nvt(dt, q, p, invM, ham, py_funct, params, 0, therm);


}// Verlet1





}// namespace libdyn
}// liblibra


