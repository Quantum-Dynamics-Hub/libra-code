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

#ifndef ANNIHILATE_H
#define ANNIHILATE_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{

void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb, MATRIX* Ra, MATRIX* Rb);
void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb);


}// namespace libcalculators

#endif // ANNIHILATE_H
