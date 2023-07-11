/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libhamiltonian_atomistic.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_ATOMISTIC_H
#define LIB_ATOMISTIC_H

#include "atomistic.h"

/// liblibra namespace
namespace liblibra{


/// libatomistic namespace
namespace libatomistic{


  void export_atomistic_objects();


}// namespace libatomistic
}// liblibra


#endif// LIB_ATOMISTIC_H
