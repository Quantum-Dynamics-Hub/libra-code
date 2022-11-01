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
  \file libnhamiltonian_generic.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_nHAMILTONIAN_GENERIC_H
#define LIB_nHAMILTONIAN_GENERIC_H


#include "nHamiltonian.h"

/// liblibra namespace
namespace liblibra{


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_generic namespace
namespace libhamiltonian_generic{


void export_nhamiltonian_generic_objects();



}// namespace libhamiltonian_generic
}// namespace libhamiltonian
}// liblibra

#endif// LIB_nHAMILTONIAN_GENERIC_H
