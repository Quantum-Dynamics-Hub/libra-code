/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libsymmetry.h
  \brief The file describes Python export function
    
*/


#ifndef LIB_SYMMETRY_H
#define LIB_SYMMETRY_H


#include "Space_Groups.h"

/// liblibra namespace
namespace liblibra{

/// libsymmetry namespace
namespace libsymmetry{


void export_symmetry_objects();


}// namespace libsymmetry
}// namespace liblibra


#endif// LIB_SYMMETRY_H
