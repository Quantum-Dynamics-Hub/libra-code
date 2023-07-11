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
  \file Potentials_elec.h
  This file defines the interfaces to the electrostatic interaction potentials
    
*/


#ifndef POTENTIALS_ELEC_H
#define POTENTIALS_ELEC_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

namespace libpot{



//------------------------- Electrostatic potentials ----------------------------------------

double Elec_Coulomb(VECTOR& ri,VECTOR& rj,     /*Inputs*/
                    VECTOR& fi,VECTOR& fj,     /*Outputs*/
                    double qi,double qj,
                    double eps,double delta);  /*Parameters*/


} //namespace libpot
}// liblibra

#endif //POTENTIALS_ELEC_H
