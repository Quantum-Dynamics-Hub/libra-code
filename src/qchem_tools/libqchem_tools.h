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
  \file libqchem_tools.h
  \brief The file describes Python export function
    
*/

#ifndef LIBQCHEM_TOOLS_H
#define LIBQCHEM_TOOLS_H

#include "Charge_Density.h"
#include "DOS.h"

/// liblibra namespace
namespace liblibra{


/// libqchem_tools namespace
namespace libqchem_tools{


void export_qchem_tools_objects();

}// namespace libqchem_tools
}// liblibra


#endif // LIBQCHEM_TOOLS_H

