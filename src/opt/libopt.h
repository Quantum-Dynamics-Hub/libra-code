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
  \file libopt.h
  \brief The file describes Python export function
    
*/

#ifndef LIBOPT_H
#define LIBOPT_H

#include "opt.h"

/// liblibra namespace
namespace liblibra{

/// libopt namespace
namespace libopt{


void export_opt_objects();


}// namespace libopt
}// liblibra

#endif // LIBOPT_H

