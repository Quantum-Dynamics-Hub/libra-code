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
/**
  \file Dynamics_Nuclear.h
  \brief The file describes the functions for nuclear (classical) dynamics
    
*/

#ifndef DYNAMICS_H
#define DYNAMICS_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../hamiltonian/libhamiltonian.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{


// Verlet.cpp
void Verlet(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params);


}// namespace libdyn
}// liblibra

#endif // DYNAMICS_H

