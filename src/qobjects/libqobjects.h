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
  \file libqobjects.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_QOBJECTS_H
#define LIB_QOBJECTS_H

#include "PrimitiveG.h"
#include "AO.h"
#include "PW.h"
#include "SD.h"

/// liblibra namespace
namespace liblibra{


/// libqobjects namespace
namespace libqobjects{


void export_qobjects_objects();


}// namespace libqobjects
}// namespace liblibra



#endif// LIB_QOBJECTS_H
