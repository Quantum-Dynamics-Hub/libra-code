/*********************************************************************************
* Copyright (C) 2018 Brendan A. Smith, Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Models_1_state.h
  \brief The file describes the 1 electronic state model potentials. They can depend on multiple 
  nuclear DOFs though
    
*/

#ifndef MODELS_1_STATE_H
#define MODELS_1_STATE_H

#include "../../math_linalg/liblinalg.h"

/// liblibra namespace
namespace liblibra{

using namespace liblinalg;


/// libhamiltonian namespace
namespace libhamiltonian{

/// libhamiltonian_model namespace
namespace libhamiltonian_model{


void model_anharmonic_1S_1D(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                          vector<double>& q, vector<double>& params);

void model_double_well_1S_1D(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia,
                             vector<CMATRIX>& dc1_dia, vector<double>& q, vector<double>& params);

void model_harmonic_1S_1D(CMATRIX& Hdia, CMATRIX& Sdia, vector<CMATRIX>& d1ham_dia, vector<CMATRIX>& dc1_dia,
                          vector<double>& q, vector<double>& params);



}// namespace libhamiltonian_model
}// namespace libhamiltonian
}// liblibra

#endif // MODEL_1_STATE_H
