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

#ifndef POTENTIALS_OOP_H
#define POTENTIALS_OOP_H

#include "../math_linalg/liblinalg.h"
#include "../math_specialfunctions/libspecialfunctions.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;
using namespace libspecialfunctions;

namespace libpot{

//------------------ Out-of-plane (OOP) potentials -------------------------------

double OOP_Fourier(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,        /*Inputs*/
                   VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,        /*Outputs*/
                   double Kijkl,double C0,double C1,double C2,int opt  /*Parameters*/
                   );

double OOP_Wilson(VECTOR& r1,VECTOR& r2,VECTOR& r3,VECTOR& r4,   /*Inputs*/
                  VECTOR& f1,VECTOR& f2,VECTOR& f3,VECTOR& f4,   /*Outputs*/
                  double Kijkl,double xi_0                       /*Parameters*/
                   );

double OOP_Harmonic(VECTOR& r0,VECTOR& r1,VECTOR& r2,VECTOR& r3,  /*Inputs*/
                    VECTOR& f0,VECTOR& f1,VECTOR& f2,VECTOR& f3,  /*Outputs*/
                    double Kijkl                                  /*Parameters*/
                   );


}// namespace libpot
}// liblibra

#endif //POTENTIALS_OOP_H
