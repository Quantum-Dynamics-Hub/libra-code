/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libnuclear.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_NUCLEAR_H
#define LIB_NUCLEAR_H


#include "Nuclear.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

/// libnuclear namespace
namespace libnuclear{


void export_Nuclear_objects();

}// namespace libdyn
}// namespace libnuclear
}// liblibra


#endif// LIB_NUCLEAR_H
