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
  \file Verlet.cpp
  \brief The file implements the velocity Verlet algorithm for nuclear dynamics integration
    
*/

#include "Dynamics.h"


/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;


/// libdyn namespace 
namespace libdyn{



void Verlet0(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params){
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
  int nadi = ham.nadi;
  int dof;
  int ham_rep = 0; // default -- assume the Hamiltonian is first computed in the diabatic representation
                   // and then will be transformed to the adiabatic in this function. 

  int act_state = 0; // default -- assume the active adiabatic surface is the ground state


  std::string key;
  boost::python::dict d = (boost::python::dict)params;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);
    if(key=="ham_rep") { ham_rep = extract<int>(d.values()[i]);   }
    if(key=="act_state"){ act_state = extract<int>(d.values()[i]); }
  }


  CMATRIX Cadi(nadi,1); Cadi.set(act_state,0, 1.0,0.0);
  
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

  p = p + ham.forces_adi(Cadi).real() * 0.5*dt;

}// Verlet0



void Verlet1(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, int entanglement_opt){
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
  int nadi = ham.nadi;
  int traj, dof;


  int ham_rep = 0; // default -- assume the Hamiltonian is first computed in the diabatic representation
                   // and then will be transformed to the adiabatic in this function. 

  int act_state = 0; // default -- assume the active adiabatic surface is the ground state



  //============= Extract optional parameters: needed for some execution scenarios =============
  double ETHD3_alpha = 1.0;
  double ETHD3_beta = 1.0;
  std::string key;
  boost::python::dict d = (boost::python::dict)params;
  for(int i=0;i<len(d.values());i++){
    key = extract<std::string>(d.keys()[i]);
    if(key=="ETHD3_alpha") { ETHD3_alpha = extract<double>(d.values()[i]);   }
    if(key=="ETHD3_beta") { ETHD3_beta = extract<double>(d.values()[i]);   }
    if(key=="ham_rep") { ham_rep = extract<int>(d.values()[i]);   }
    if(key=="act_state"){ act_state = extract<int>(d.values()[i]); }
  }


//  vector<int> t1(ndof, 0); for(int i=0;i<ndof;i++){  t1[i] = i; }
//  vector<int> t2(1,0);
//  vector<int> t3(2,0);

//  CMATRIX Cadi(nadi,1); Cadi.set(act_state,0, 1.0,0.0);
  MATRIX gamma(ndof, ntraj);
//  MATRIX f(ndof, 1);

  

  p = p + ham.forces_adi(act_state).real() * 0.5 * dt;


  if(entanglement_opt==22){
    gamma = ETHD3_friction(q, p, invM, ETHD3_alpha, ETHD3_beta);
  }


  // For efficiency, switch to the element-wise multiplication
  for(traj=0; traj<ntraj; traj++){
    for(dof=0; dof<ndof; dof++){  

      q.add(dof, traj,  invM.get(dof,0) * p.get(dof,traj) * dt ); 

      if(entanglement_opt==22){
        q.add(dof, traj,  invM.get(dof,0) * gamma.get(dof,traj) * dt ); 
      }

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
  else if(entanglement_opt==22){  ham.add_ethd3_adi(q, p, invM, ETHD3_alpha, ETHD3_beta, 1);  }
  else{
    cout<<"ERROR in Verlet1: The entanglement option = "<<entanglement_opt<<" is not avaialable\n";
    exit(0);
  }

  p = p + ham.forces_adi(act_state).real() * 0.5 * dt;




}// Verlet1


void Verlet1(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params){
/**
  \brief Regular velocity Verlet for an ensemble of (uncoupled trajectories)
*/

  Verlet1(dt, q, p, invM, ham, py_funct, params, 0);


}// Verlet1





}// namespace libdyn
}// liblibra


