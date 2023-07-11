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
  \file libtimer.h
  \brief The file describes Python export function
    
*/

#ifndef LIBTIMER_H
#define LIBTIMER_H

#include "Timer.h"

/// liblibra namespace
namespace liblibra{

void export_timer_objects();


}// namespace libtimer

#endif // LIBTIMER_H

