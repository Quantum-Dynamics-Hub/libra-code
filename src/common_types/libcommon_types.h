/*********************************************************************************
* Copyright (C) 2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libcommon_types.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_COMMON_TYPES_H
#define LIB_COMMON_TYPES_H

#include "Excitation.h"

/// liblibra namespace
namespace liblibra{

/// libcommon_types namespace
namespace libcommon_types{


void export_common_types_objects();


}// namespace libcommon_types
}// namespace liblibra


#endif// LIB_COMMON_TYPES_H
