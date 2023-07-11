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


void update_Hamiltonian_variables(dyn_control_params& prms, dyn_variables& dyn_var, 
                                  nHamiltonian& ham, nHamiltonian& ham_prev,
                                  bp::object py_funct, bp::object model_params, int update_type);
void update_Hamiltonian_variables(bp::dict prms, dyn_variables& dyn_var, 
                                  nHamiltonian& ham, nHamiltonian& ham_prev,
                                  bp::object py_funct, bp::object model_params, int update_type);



}// namespace libdyn
}// liblibra

#endif // DYN_HAM_H

