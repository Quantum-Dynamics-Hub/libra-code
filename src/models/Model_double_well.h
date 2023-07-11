/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Model_SAC.h
  \brief The file describes the functions for computing SAC (single avoided crossing) Hamiltonian and its derivatives
    
*/

#ifndef MODEL_DOUBLE_WELL_H
#define MODEL_DOUBLE_WELL_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libmodels namespace
namespace libmodels{


void double_well_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params_);
boost::python::list double_well_Ham(double x, boost::python::list params_);


}// namespace libmodels
}// liblibra

#endif // MODEL_DOUBLE_WELL_H
