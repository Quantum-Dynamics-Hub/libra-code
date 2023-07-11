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
  \file libuniverse.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_UNIVERSE_H
#define LIB_UNIVERSE_H

#include "Element.h"
#include "Universe.h"

/// liblibra namespace
namespace liblibra{

/// libchemobjects namespace
namespace libchemobjects{

/// libuniverse namespace
namespace libuniverse{

void export_Universe_objects();


}// namespace libuniverse
}// namespace libchemobjects
}// liblibra

#endif // LIB_UNIVERSE_H
