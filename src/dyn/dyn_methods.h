/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_methods.h
  \brief The file describes the stand-alone methods for dynamics
    
*/

#ifndef DYN_METHODS_H
#define DYN_METHODS_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{



///================  In dyn_methods_dish.cpp  ===================================
/*
int ida(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi);
int dish(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib,
          int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2);
int dish(Electronic& el, Nuclear& mol, Hamiltonian& ham, 
          MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2);
*/




}// namespace libdyn
}// liblibra

#endif // DYN_METHODS_H

