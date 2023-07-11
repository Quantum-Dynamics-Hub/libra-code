/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/

#ifndef MODEL_SIN_2D_H
#define MODEL_SIN_2D_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


namespace libmodels{

void sin_2D_Ham(double x, double y, MATRIX* H, 
                MATRIX* dH1,  MATRIX* dH2, 
                MATRIX* d2H1, MATRIX* d2H2, vector<double>& params);
boost::python::list sin_2D_Ham(double x, double y, boost::python::list params_);


}// namespace libhamiltonian
}// liblibra

#endif // MODEL_SIN_2D_H
