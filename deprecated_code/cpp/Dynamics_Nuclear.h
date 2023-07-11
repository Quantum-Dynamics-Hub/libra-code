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

#ifndef DYNAMICS_NUCLEAR_H
#define DYNAMICS_NUCLEAR_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../hamiltonian/libhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
#include "thermostat/libthermostat.h"

/// liblibra namespace
namespace liblibra{

using namespace libhamiltonian;

/// libdyn namespace
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace libthermostat;


// Dynamics_Nuclear.cpp
void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,int opt);
void propagate_nuclear(double dt,Nuclear* mol,Electronic* el,Hamiltonian* ham,Thermostat* therm, int opt);


}// namespace libdyn
}// liblibra

#endif // DYNAMICS_NUCLEAR_H

