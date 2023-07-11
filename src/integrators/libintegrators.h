/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libintegrators.h
  \brief The file describes Python export function
    
*/

#ifndef LIBINTEGRATORS_H
#define LIBINTEGRATORS_H

#include "integrators.h"


/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libintegrators{


void export_integrators_objects();


}// namespace libintegrators
}// liblibra

#endif // LIBINTEGRATORS_H

