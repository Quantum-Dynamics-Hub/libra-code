/*********************************************************************************
* Copyright (C) 2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file NPI.h
  \brief The file contains prototypes for the norm-preserving interpolation (NPI) method to compute NACs
    
*/

#ifndef NPI_H
#define NPI_H


#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libcalculators namespace
namespace libcalculators{

MATRIX nac_npi(MATRIX& St, double dt);


}// namespace libcalculators
}// liblibra

#endif // NPI_H
