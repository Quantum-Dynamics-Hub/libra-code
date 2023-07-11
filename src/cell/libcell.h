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
  \file libcell.h
  \brief The file describes Python export function
    
*/

#ifndef LIBCELL_H
#define LIBCELL_H

#include "Cell.h"
#include "NList.h"

/// liblibra namespace
namespace liblibra{

/// libcell namespace
namespace libcell{

void export_Cell_objects();

}// namespace libcell
}// liblibra

#endif //
