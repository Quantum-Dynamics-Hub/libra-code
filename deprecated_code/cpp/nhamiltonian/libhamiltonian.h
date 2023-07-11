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
  \file libhamiltonian.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_HAMILTONIAN_H
#define LIB_HAMILTONIAN_H

#include "nHamiltonian_Generic/libnhamiltonian_generic.h"
//#include "Hamiltonian_Generic/libhamiltonian_generic.h"
//#include "Hamiltonian_Model/libhamiltonian_model.h"
//#include "Hamiltonian_Atomistic/libhamiltonian_atomistic.h"
//#include "Hamiltonian_Extern/libhamiltonian_extern.h"

/// liblibra namespace
namespace liblibra{


namespace libnhamiltonian{

void export_Hamiltonian_objects();


}// namespace libnhamiltonian
}// liblibra

#endif// LIB_HAMILTONIAN_H
