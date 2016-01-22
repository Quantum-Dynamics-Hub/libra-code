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
  \file libio.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_IO_H
#define LIB_IO_H


#include "io.h"

/// libio namespace
namespace libio{


void export_io_objects();


}// namespace libio


#endif// LIB_IO_H
