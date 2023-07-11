/*********************************************************************************
* Copyright (C) 2015-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#include "VECTOR.h"
#include "QUATERNION.h"
#include "MATRIX3x3.h"
#include "MATRIX.h"

/// liblibra namespace
namespace liblibra{

/// liblinalg namespace
namespace liblinalg{



void MATRIX_TO_QUATERNION(MATRIX&,QUATERNION&);
void QUATERNION_TO_MATRIX(QUATERNION&,MATRIX&);
void MATRIX_TO_QUATERNION(MATRIX3x3&,QUATERNION&);
void QUATERNION_TO_MATRIX(QUATERNION&,MATRIX3x3&);


}// namespace liblinalg
}// namespace liblibra

