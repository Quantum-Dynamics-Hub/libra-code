/*********************************************************************************
* Copyright (C) 2015-2022 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Model_DAC.h
  \brief The file describes the functions for computing DAC (Double Avoided Crossing) Hamiltonian and its derivatives
    
*/

#ifndef MODEL_DAC_H
#define MODEL_DAC_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;

/// libmodels namespace
namespace libmodels{


void model_DAC(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
               vector<double> &q, vector<double>& params);

void DAC_Ham(double x, MATRIX* H, MATRIX* dH, MATRIX* d2H, vector<double>& params_);
boost::python::list DAC_Ham(double x, boost::python::list params_);


}// namespace libmodels
}// liblibra

#endif // MODEL_DAC_H
