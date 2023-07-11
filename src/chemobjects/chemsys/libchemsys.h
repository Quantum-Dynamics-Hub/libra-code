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
  \file libchemsys.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_CHEMSYS_H
#define LIB_CHEMSYS_H


#include "System.h"

/// liblibra namespace
namespace liblibra{


/// libchemobjects namespace
namespace libchemobjects{

/// libchemsys namespace
namespace libchemsys{


void export_Chemsys_objects();


}// namespace libchemsys
}// namespace libchemobjects

}// liblibra


#endif// LIB_CHEMSYS_H
