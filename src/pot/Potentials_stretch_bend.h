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

#ifndef POTENTIALS_STRETCH_BEND_H
#define POTENTIALS_STRETCH_BEND_H


#include "../math_linalg/liblinalg.h"
#include "../Units.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{



//------------------- Stretch-bend potentials --------------------------------
double Stretch_Bend_Harmonic(VECTOR& r1,VECTOR& r2,VECTOR& r3,         /* Inputs */
                             VECTOR& f1,VECTOR& f2,VECTOR& f3,         /* Outputs*/
                             double k_ijk,double k_kji, double theta_0,
                             double r_ij0,double r_kj0                 /* Parameters*/
                            );


}// namespace libpot
}// liblibra

#endif //POTENTIALS_STRETCH_BEND_H
