/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef MODEL_MARCUS_H
#define MODEL_MARCUS_H

#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libhamiltonian{
namespace libhamiltonian_model{


void Marcus_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params_);
boost::python::list Marcus_Ham(double x, boost::python::list params_);


}// namespace libhamiltonian_model
}// namespace libhamiltonian
}// liblibra

#endif // MODEL_MARCUS
