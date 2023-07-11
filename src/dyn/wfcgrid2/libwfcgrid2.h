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
  \file libwfcgrid2.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_WFCGRID2_H
#define LIB_WFCGRID2_H


#include "Wfcgrid2.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid2{


void export_Wfcgrid2_objects();


}// namespace libwfcgrid2
}// namespace libdyn
}// liblibra



#endif// LIB_WFCGRID2_H
