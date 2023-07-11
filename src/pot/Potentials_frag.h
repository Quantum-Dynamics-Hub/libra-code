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

#ifndef POTENTIALS_FRAG_H
#define POTENTIALS_FRAG_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

namespace libpot{


//------------------------- Fragment-Fragment potentials ------------------------------------
double Gay_Berne(VECTOR& ri,VECTOR& rj,VECTOR& ui,VECTOR& uj,          /*Inputs*/
                 VECTOR& fi,VECTOR& fj,VECTOR& ti,VECTOR& tj,          /*Outputs*/
                 double di, double dj,double li,double lj,
                 double e0,double rat,double dw,double mu,double nu);  /*Parameters*/

double Girifalco12_6(VECTOR& ri,VECTOR& rj,
                     VECTOR& fi,VECTOR& fj,
                     double a,double alp,double bet
                    );

}//namespace libpot
}// liblibra

#endif // POTENTIALS_FRAG_H
