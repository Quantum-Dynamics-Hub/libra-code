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
  \file libcontrol_parameters.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_CONTROL_PARAMETERS_H
#define LIB_CONTROL_PARAMETERS_H

#include "Control_Parameters.h"

/// liblibra namespace
namespace liblibra{


/// libcontrol_parameters namespace
namespace libcontrol_parameters{


void export_Control_Parameters_objects();


}// namespace libcontrol_parameters
}// namespace liblibra


#endif// LIB_CONTROL_PARAMETERS_H
