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
  \file integrators.h
  \brief 
    
*/

#ifndef INTEGRATORS_H
#define INTEGRATORS_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;
//using namespace boost::python;

/// libdyn namespace
namespace libintegrators{

namespace bp = boost::python;

CMATRIX RK4(CMATRIX& q, double dt, bp::object compute_derivatives, bp::object function_params);


}// namespace libintegrators

}// liblibra

#endif // INTEGRATORS_H

