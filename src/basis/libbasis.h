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
  \file libbasis.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_BASIS_H
#define LIB_BASIS_H

#include "Basis.h"
#include "basisset.h"

/// liblibra namespace
namespace liblibra{



/// libbasis namespace
namespace libbasis{


void export_basis_objects();


}// namespace libbasis
}// namespace liblibra



#endif// LIB_BASIS_H
