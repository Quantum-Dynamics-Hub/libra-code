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
  \file liboperators.h
  \brief The file describes Python export function
    
*/

#ifndef LIBOPERATORS_H
#define LIBOPERATORS_H


#include "Operators.h"

/// libmmath namespace
namespace libmmath{

/// liboperators namespace
namespace liboperators{


void export_Operators_objects();


}// namespace liboperators
}// namespace libmmath

#endif// LIBOPERATORS_H
