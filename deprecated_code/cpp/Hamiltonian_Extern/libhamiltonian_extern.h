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
  \file libhamiltonian_extern.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_HAMILTONIAN_EXTERN_H
#define LIB_HAMILTONIAN_EXTERN_H

#include "Hamiltonian_Extern.h"

/// liblibra namespace
namespace liblibra{


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_extern namespace
namespace libhamiltonian_extern{

void export_hamiltonian_extern_objects();


}// namespace libhamiltonian_extern
}// namespace libhamiltonian
}// liblibra

#endif// LIB_HAMILTONIAN_EXTERN_H
