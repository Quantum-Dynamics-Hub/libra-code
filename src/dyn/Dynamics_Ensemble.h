/*********************************************************************************
* Copyright (C) 2015 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Dynamics_Ensemble.h
  \brief The file describes the functions for dynamics of ensemble of trajectories
    
*/

#ifndef DYNAMICS_ENSEMBLE_H
#define DYNAMICS_ENSEMBLE_H

// External dependencies
#include "../mmath/libmmath.h"
#include "../hamiltonian/libhamiltonian.h"

// Dynamics classes
#include "nuclear/libnuclear.h"
#include "electronic/libelectronic.h"
#include "thermostat/libthermostat.h"
#include "ensemble/libensemble.h"

#include "Dynamics_Nuclear.h"


using namespace libhamiltonian;

/// libdyn namespace 
namespace libdyn{

using namespace libnuclear;
using namespace libelectronic;
using namespace libthermostat;
using namespace libensemble;



void propagate_ensemble(double dt,Ensemble* ens,int opt);
void propagate_ensemble(double dt,Ensemble& ens,int opt);


}// namespace libdyn

#endif // DYNAMICS_ENSEMBLE_H
