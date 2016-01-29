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
  \file libcontrol_parameters.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_CONTROL_PARAMETERS_H
#define LIB_CONTROL_PARAMETERS_H

#include "Control_Parameters.h"

/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_atomistic namespace
namespace libhamiltonian_atomistic{

/// libhamiltonian_qm namespace
namespace libhamiltonian_qm{

/// libcontrol_parameters namespace
namespace libcontrol_parameters{


void export_Control_Parameters_objects();


}// namespace libcontrol_parameters
}// namespace libhamiltonian_qm
}// namespace libhamiltonian_atomistic
}// namespace libhamiltonian


#endif// LIB_CONTROL_PARAMETERS_H
