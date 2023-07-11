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

#ifndef SWITCHING_FUNCTIONS_H
#define SWITCHING_FUNCTIONS_H


#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;



namespace libpot{

//------------------ Switching functions --------------------------
void SWITCH(VECTOR& r1,VECTOR&r2, double R_on,double R_off,double& SW,VECTOR& dSW);
boost::python::list SWITCH(VECTOR r1,VECTOR r2, double R_on,double R_off);

void DOUBLE_SWITCH(double x,double a,double eps,double& SW,double& dSW);
boost::python::list DOUBLE_SWITCH(double x,double a,double eps);



}// namespace libpot
}// liblibra


#endif // SWITCHING_FUNCTIONS_H
