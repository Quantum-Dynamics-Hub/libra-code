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
  \file Energy_Nuclear.h
  \brief The file defines functions for nuclear interaction energy/force calculations
    
*/

#ifndef ENERGY_NUCLEAR_H
#define ENERGY_NUCLEAR_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libcalculators namespace
namespace libcalculators{

double energy_nucl(vector<VECTOR>& R, vector<double>& Zeff);
double energy_nucl(vector<VECTOR>& R, vector<double>& Zeff, vector<VECTOR>& G);


}// namespace libcalculators

}// liblibra

#endif // ENERGY_NUCLEAR_H
