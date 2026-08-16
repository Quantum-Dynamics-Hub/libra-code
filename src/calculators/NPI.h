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

/**
 * Compute the interval-averaged real time-derivative coupling using the
 * Meek-Levine norm-preserving interpolation.
 *
 * St must be a finite, square, phase-matched orthogonal overlap matrix with
 * positive determinant, and dt must be finite and positive. Invalid inputs
 * raise std::invalid_argument with guidance for correcting the overlap.
 */
MATRIX nac_npi(MATRIX& St, double dt);


}// namespace libcalculators
}// liblibra

#endif // NPI_H
