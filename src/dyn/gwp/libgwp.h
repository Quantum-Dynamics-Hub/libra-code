/*********************************************************************************
* Copyright (C) 2015-2020 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libgwp.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_GWP_H
#define LIB_GWP_H


#include "gwp.h"

/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libdyn{

/// libgwp namespace
namespace libgwp{


void export_gwp_objects();


}// namespace libgwp
}// namespace libdyn
}// liblibra



#endif// LIB_GWP_H
