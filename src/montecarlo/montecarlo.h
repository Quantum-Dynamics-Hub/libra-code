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
  \file montecarlo.h
  \brief The file describes Monte Carlo sampling procedures
    
*/

#ifndef MONTECARLO_H
#define MONTECARLO_H


// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../math_random/librandom.h"
#include "../io/libio.h"



/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace librandom;
using namespace libio;

/// libmontecarlo namespace 
namespace libmontecarlo{

namespace bp = boost::python;

vector<MATRIX> metropolis_gau(Random& rnd, bp::object target_distribution, MATRIX& dof, bp::object distribution_params, 
                              int sample_size, int start_sampling, double gau_var);



}// namespace libmontecarlo
}// liblibra

#endif // MONTECARLO_H
