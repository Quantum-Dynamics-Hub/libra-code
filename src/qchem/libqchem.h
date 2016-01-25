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
  \file libqchem.h
  \brief The file describes Python export function
    
*/

#ifndef LIBQCHEM_H
#define LIBQCHEM_H

#include "molint/libmolint.h"
#include "qobjects/libqobjects.h"
#include "basis/libbasis.h"

/// libqchem namespace
namespace libqchem{


void export_Qchem_objects();



}// namespace libqchem

#endif // LIBQCHEM_H

