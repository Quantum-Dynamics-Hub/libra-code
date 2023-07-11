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
  \file libthermostat.h
  \brief The file describes Python export function
    
*/

#ifndef LIB_THERMOSTAT_H
#define LIB_THERMOSTAT_H


#include "Thermostat.h"

/// liblibra namespace
namespace liblibra{

/// libdyn namespace
namespace libdyn{

/// libthermostat namespace
namespace libthermostat{


void export_Thermostat_objects();


}// namespace libdyn
}// namespace libthermostat
}// liblibra


#endif// LIB_THERMOSTAT_H
