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
  \file libmolint.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_MOLINT_H
#define LIB_MOLINT_H


#include "A_coefficients.h"
#include "Overlaps.h"
#include "Moments.h"
#include "Pseudopotential.h"
#include "Multipoles.h"

#include "Integral_Kinetic.h"
#include "Integral_Nuclear_Attraction.h"
#include "Integral_Electron_Repulsion.h"
#include "Integral_Derivative_Couplings.h"
#include "Integral_Approx1.h"


/// liblibra namespace
namespace liblibra{


/// libmolint namespace
namespace libmolint{


void export_molint_objects();


}// namespace libmolint
}// namespace liblibra



#endif// LIB_MOLINT_H
