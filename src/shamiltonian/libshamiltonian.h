/*********************************************************************************
* Copyright (C) 2025 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libshamiltonian_generic.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_sHAMILTONIAN_H
#define LIB_sHAMILTONIAN_H


#include "sHamiltonian.h"

/// liblibra namespace
namespace liblibra{


/// libshamiltonian namespace
namespace libshamiltonian{


void export_shamiltonian_objects();



}// namespace libshamiltonian
}// liblibra

#endif// LIB_sHAMILTONIAN_H
