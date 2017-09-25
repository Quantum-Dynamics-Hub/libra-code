/*********************************************************************************
* Copyright (C) 2015-2017 Brendan A. Smith , Alexey V. Akimov 
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef MODEL_DOUBLE_WELL_H
#define MODEL_DOUBLE_WELL_H

#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_model namespace
namespace libhamiltonian_model{


void double_well_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params);
boost::python::list double_well_Ham(double x, boost::python::list params_);


}// namespace libhamiltonian_model
}// namespace libhamiltonian
}// liblibra

#endif // MODEL_DOUBLE_WELL_H
