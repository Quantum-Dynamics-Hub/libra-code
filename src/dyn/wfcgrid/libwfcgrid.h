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
  \file libwfcgrid.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_WFCGRID_H
#define LIB_WFCGRID_H


#include "Wfcgrid.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

/// libwfcgrid namespace
namespace libwfcgrid{


void export_Wfcgrid_objects();


}// namespace libwfcgrid
}// namespace libdyn
}// liblibra



#endif// LIB_WFCGRID_H
