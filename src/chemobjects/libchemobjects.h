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
  \file libchemobjects.h
  \brief The file describes Python export function
    
*/

#ifndef LIBCHEMOBJECTS_H
#define LIBCHEMOBJECTS_H

#include "universe/libuniverse.h"
#include "mol/libmol.h"
#include "chemsys/libchemsys.h"

/// liblibra namespace
namespace liblibra{


/// libchemobjects namespace
namespace libchemobjects{


void export_chemobjects_objects();


}// namespace libchemobjects
}// liblibra

#endif // LIBCHEMOBJECTS_H

