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
  \file libbarostat.h
  \brief The file describes Python export function
    
*/


#ifndef LIB_BAROSTAT_H
#define LIB_BAROSTAT_H


#include "Barostat.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

/// libbarostat namespace
namespace libbarostat{


void export_Barostat_objects();


}// namespace libdyn
}// namespace libbarostat
}// liblibra


#endif// LIB_BAROSTAT_H
