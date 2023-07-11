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
  \file opt.h
  \brief The file describes various optimization procedures
    
*/

#ifndef OPT_H
#define OPT_H


// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"



/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libio;

/// libopt namespace 
namespace libopt{

namespace bp = boost::python;

MATRIX grad_descent(bp::object grad_function, MATRIX& dof, bp::object funct_params, double grad_tol, double step_size, int max_steps);


}// namespace libopt
}// liblibra

#endif // OPT_H
