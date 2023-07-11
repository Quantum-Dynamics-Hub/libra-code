/*********************************************************************************
* Copyright (C) 2022 Matthew Dutra, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libqtag.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_QTAG_H
#define LIB_QTAG_H


#include "qtag.h"

/// liblibra namespace
namespace liblibra{


/// libdyn namespace
namespace libdyn{

/// libqtag namespace
namespace libqtag{


void export_qtag_objects();


}// namespace libqtag
}// namespace libdyn
}// liblibra



#endif// LIB_QTAG_H
