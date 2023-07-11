/*********************************************************************************
* Copyright (C) 2016-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libutil.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_UTIL_H
#define LIB_UTIL_H

#include "util.h"

/// liblibra namespace
namespace liblibra{


/// libutil namespace
namespace libutil{


void export_util_objects();


}// namespace libutil
}// liblibra


#endif// LIB_UTIL_H
