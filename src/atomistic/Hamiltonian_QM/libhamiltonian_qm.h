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
  \file libhamiltonian_qm.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_HAMILTONIAN_QM_H
#define LIB_HAMILTONIAN_QM_H


#include "Hamiltonian_QM.h"
#include "SCF.h"

/// liblibra namespace
namespace liblibra{

namespace libatomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{


void export_hamiltonian_qm_objects();



}// namespace libhamiltonian_qm
}// namespace libatomistic
}// liblibra

#endif// LIB_HAMILTONIAN_QM_H
