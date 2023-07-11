/*********************************************************************************
* Copyright (C) 2019-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libheom.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_HEOM_H
#define LIB_HEOM_H


#include "heom.h"

namespace liblibra{
namespace libdyn{
namespace libheom{




void export_heom_objects();


}// namespace libheom
}// namespace libdyn
}// liblibra



#endif// LIB_HEOM_H
