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

#ifndef POTENTIALS_VDW_H
#define POTENTIALS_VDW_H

#include "Switching_functions.h"
#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libpot{


//--------------------------- Vdw potentials -------------------------------------------

double Vdw_LJ(VECTOR& ri,VECTOR& rj,          /*Inputs*/
              VECTOR& fi,VECTOR& fj,          /*Outputs*/
              double sigma, double espilon);  /*Parameters*/

double Vdw_Buffered14_7(VECTOR& ri,VECTOR& rj,          /*Inputs*/
                        VECTOR& fi,VECTOR& fj,          /*Outputs*/
                        double sigma, double espilon);  /*Parameters*/

double Vdw_Morse(VECTOR& ri,VECTOR& rj,            /*Inputs*/
                 VECTOR& fi,VECTOR& fj,            /*Outputs*/
                 double D, double r0,double alp);  /*Parameters*/


} // namespace libpot
}//liblibra

#endif //POTENTIALS_VDW_H
