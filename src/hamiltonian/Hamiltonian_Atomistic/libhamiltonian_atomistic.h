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
  \file libhamiltonian_atomistic.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_HAMILTONIAN_ATOMISTIC_H
#define LIB_HAMILTONIAN_ATOMISTIC_H

#include "Hamiltonian_Atomistic.h"

/// liblibra namespace
namespace liblibra{


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_atomistic namespace
namespace libhamiltonian_atomistic{


void export_hamiltonian_atomistic_objects();


}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian
}// liblibra

#endif// LIB_HAMILTONIAN_ATOMISTIC_H
