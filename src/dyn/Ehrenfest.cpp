/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Ehrenfest.cpp
  \brief The file implements all about Ehrenfest calculations
    
*/

#include "electronic/libelectronic.h"
#include "Dynamics.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace 
namespace libdyn{

using namespace libelectronic;


void Ehrenfest0(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, nHamiltonian& ham, bp::object py_funct, bp::object params, int rep){
/**
  \brief One step of the Ehrenfest algorithm for electron-nuclear DOFs for one trajectory

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] nadi x nadi or ndia x ndia matrix containing the electronic coordinates
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) - its internal variables
  (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] rep The representation to run the calculations: 0 - diabatic, 1 - adiabatic

*/

 
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



void Ehrenfest1(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, nHamiltonian& ham, bp::object py_funct, bp::object params, int rep){
/**
  \brief One step of the Ehrenfest algorithm for electron-nuclear DOFs for an ensemble of trajectories

  \param[in] Integration time step
  \param[in,out] q [Ndof x Ntraj] nuclear coordinates. Change during the integration.
  \param[in,out] p [Ndof x Ntraj] nuclear momenta. Change during the integration.
  \param[in] invM [Ndof  x 1] inverse nuclear DOF masses. 
  \param[in,out] nadi x nadi or ndia x ndia matrix containing the electronic coordinates
  \param[in] ham Is the Hamiltonian object that works as a functor (takes care of all calculations of given type) - its internal variables
  (well, actually the variables it points to) are changed during the compuations
  \param[in] py_funct Python function object that is called when this algorithm is executed. The called Python function does the necessary 
  computations to update the diabatic Hamiltonian matrix (and derivatives), stored externally.
  \param[in] params The Python object containing any necessary parameters passed to the "py_funct" function when it is executed.
  \param[in] rep The representation to run the calculations: 0 - diabatic, 1 - adiabatic

*/

  cout<<"ERROR: the implementation is in process. Exiting...\n";
  exit(0);
 
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







}// namespace libdyn
}// liblibra
