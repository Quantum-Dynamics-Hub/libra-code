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
  \file libhamiltonian_mm.h
  \brief The file describes Python export function
    
*/


#ifndef LIB_HAMILTONIAN_MM_H
#define LIB_HAMILTONIAN_MM_H

#include "Hamiltonian_MM.h"


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_atomistic namespace
namespace libhamiltonian_atomistic{

/// libhamiltonian_mm namespace
namespace libhamiltonian_mm{


void export_Hamiltonian_MM_objects();



}// namespace libhamiltonian_mm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


#endif// LIB_HAMILTONIAN_MM_H
