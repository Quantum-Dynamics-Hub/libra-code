/*********************************************************************************
* Copyright (C) 2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_methods.h
  \brief The file describes the stand-alone methods for dynamics
    
*/

#ifndef DYN_METHODS_H
#define DYN_METHODS_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"
#include "dyn_variables.h"
#include "../nhamiltonian/libnhamiltonian.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libnhamiltonian;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{



///================  In dyn_methods_dish.cpp  ===================================

vector<int> decoherence_event(MATRIX& coherence_time, MATRIX& coherence_interval, int decoherence_event_option, Random& rnd);
vector<int> decoherence_event(MATRIX& coherence_time, MATRIX& coherence_interval, Random& rnd);

vector<int> dish_hop_proposal(vector<int>& act_states, CMATRIX& Coeff, 
  MATRIX& coherence_time, vector<MATRIX>& decoherence_rates, Random& rnd);

void dish_project_out_collapse(vector<int>& old_states, vector<int>& proposed_states, vector<int>& new_states, 
  CMATRIX& Coeff, MATRIX& coherence_time, int collapse_option);


//================= In dyn_methods_qtsh.cpp =======================================

MATRIX compute_dkinemat(dyn_variables& dyn_var, nHamiltonian& ham);




}// namespace libdyn
}// liblibra

#endif // DYN_METHODS_H

