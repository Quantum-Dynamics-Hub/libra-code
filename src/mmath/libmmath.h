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

#ifndef LIBMMATH_H
#define LIBMMATH_H

#include "linalg/liblinalg.h"
#include "specialfunctions/libspecialfunctions.h"
#include "graph/libgraph.h"
#include "operators/liboperators.h"
#include "random/librandom.h"
#include "data/libdata.h"
#include "ann/libann.h"
#include "meigen/libmeigen.h"
#include "symmetry/libsymmetry.h"

#include "Timer.h"


namespace libmmath{

void export_Mathematics_objects();


}// namespace libmmath

#endif // LIBMMATH_H

