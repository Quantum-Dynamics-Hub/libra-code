/*********************************************************************************
* Copyright (C) 2018 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file libmontecarlo.h
  \brief The file describes Python export function
    
*/

#ifndef LIBMONTECARLO_H
#define LIBMONTECARLO_H

#include "montecarlo.h"

/// liblibra namespace
namespace liblibra{

/// libmontecarlo namespace
namespace libmontecarlo{


void export_montecarlo_objects();


}// namespace libmontecarlo
}// liblibra

#endif // LIBMONTECARLO_H

