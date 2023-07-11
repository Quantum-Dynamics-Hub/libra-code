/*********************************************************************************
* Copyright (C) 2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_ham.h
  \brief The file describes the functions for updating Hamiltonian variables in response 
  to changed dynamical variables
    
*/

#ifndef DYN_HAM_H
#define DYN_HAM_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../nhamiltonian/libnhamiltonian.h"
#include "../calculators/libcalculators.h"
#include "../io/libio.h"
#include "dyn_control_params.h"
#include "dyn_variables.h"

/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libnhamiltonian;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{


void update_Hamiltonian_variables(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, 
                                  bp::object py_funct, bp::object model_params, int update_type);
void update_Hamiltonian_variables(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, 
                                  bp::object py_funct, bp::object model_params, int update_type);

void update_Hamiltonian_q(dyn_control_params& prms, MATRIX& q, nHamiltonian& ham, bp::object py_funct, bp::object model_params);
void update_Hamiltonian_q(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params);
void update_Hamiltonian_q(bp::dict prms, MATRIX& q, nHamiltonian& ham, bp::object py_funct, bp::object model_params);
void update_Hamiltonian_q(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params);


void update_Hamiltonian_q_ethd(dyn_control_params& prms, MATRIX& q, MATRIX& p, nHamiltonian& ham, 
                               bp::object py_funct, bp::object model_params, MATRIX& invM);
void update_Hamiltonian_q_ethd(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params);
void update_Hamiltonian_q_ethd(bp::dict prms, MATRIX& q, MATRIX& p, nHamiltonian& ham, 
                               bp::object py_funct, bp::object model_params, MATRIX& invM);
void update_Hamiltonian_q_ethd(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham, bp::object py_funct, bp::object model_params);


void update_Hamiltonian_p(dyn_control_params& prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM);
void update_Hamiltonian_p(dyn_control_params& prms, dyn_variables& dyn_var, nHamiltonian& ham);
void update_Hamiltonian_p(bp::dict prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM);
void update_Hamiltonian_p(bp::dict prms, dyn_variables& dyn_var, nHamiltonian& ham);


void update_nacs(dyn_control_params& prms, nHamiltonian& ham);


}// namespace libdyn
}// liblibra

#endif // DYN_HAM_H

