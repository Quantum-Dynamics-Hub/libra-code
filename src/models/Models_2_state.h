/*********************************************************************************
* Copyright (C) 2018-2022 Brendan A. Smith, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Models_2_state.h
  \brief The file describes the 2 electronic state model potentials. They can depend on multiple 
  nuclear DOFs though
    
*/

#ifndef MODELS_2_STATE_H
#define MODELS_2_STATE_H

#include "../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libmodels namespace
namespace libmodels{


vector<double> set_params_2S_1D_sin(std::string model);

void model_2S_1D_sin(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                     vector<double>& q, vector<double>& params);


vector<double> set_params_2S_2D_sin(std::string model);

void model_2S_2D_sin(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                     vector<double>& q, vector<double>& params);



vector<double> set_params_2S_1D_tanh(std::string model);

void model_2S_1D_tanh(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                      vector<double>& q, vector<double>& params);

}// namespace libmodels
}// liblibra

#endif // MODELS_2_STATE_H
