/*********************************************************************************
* Copyright (C) 2022 Matthew Dutra, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file qtag.h
  \brief This file defines the functions need for QTAG method

*/

#ifndef QTAG_H
#define QTAG_H


#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libdyn namespace
namespace libdyn{

/// libqtag namespace
namespace libqtag{


///=============== (qtag.cpp) ===================

double qtag_momentum(MATRIX& q, MATRIX& p, MATRIX& alp, MATRIX& s, CMATRIX& Coeff);




}// namespace libqtag
}// namespace libdyn
}// liblibra

#endif  // QTAG_H
