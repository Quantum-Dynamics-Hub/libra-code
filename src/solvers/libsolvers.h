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
 \file libsolvers.h
 \brief The file that exprots libsolvers objects to Python
        
*/


#ifndef LIB_SOLVERS_H
#define LIB_SOLVERS_H

#include "DIIS.h"

/// liblibra namespace
namespace liblibra{


/// libsolvers namespace
namespace libsolvers{


void export_solvers_objects();


}// namespace libsolvers
}// liblibra

#endif// LIB_SOLVERS_H
