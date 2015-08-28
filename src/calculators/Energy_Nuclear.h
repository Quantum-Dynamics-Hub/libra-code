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

#ifndef ENERGY_NUCLEAR_H
#define ENERGY_NUCLEAR_H

#include "../mmath/libmmath.h"
using namespace libmmath;

namespace libcalculators{

double energy_nucl(vector<VECTOR>& R, vector<double>& Zeff);


}// namespace libcalculators

#endif // ENERGY_NUCLEAR_H
