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
  \file libhamiltonian_mm.h
  \brief The file describes Python export function
    
*/


#ifndef LIB_HAMILTONIAN_MM_H
#define LIB_HAMILTONIAN_MM_H

#include "Hamiltonian_MM.h"
#include "Interactions.h"
#include "Interactions_2_Body.h"
#include "Interactions_3_Body.h"
#include "Interactions_4_Body.h"

/// liblibra namespace
namespace liblibra{


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
}// liblibra

#endif// LIB_HAMILTONIAN_MM_H
