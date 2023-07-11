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
  \file Annihilate.h
  \brief The file describes functions for spin annihilation - which may be needed in unrestricted (spin-polarized) calculations
    
*/

#ifndef ANNIHILATE_H
#define ANNIHILATE_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libcalculators namespace
namespace libcalculators{

void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb, MATRIX* Ra, MATRIX* Rb);
void annihilate(int Na, int Nb, MATRIX* Pa, MATRIX* Pb);


}// namespace libcalculators
}// liblibra

#endif // ANNIHILATE_H
