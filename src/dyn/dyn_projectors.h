/*********************************************************************************
* Copyright (C) 2015-2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_projectors.h
  \brief The header for dyn_projectors.cpp
    
*/

#ifndef DYN_PROJECTORS_H
#define DYN_PROJECTORS_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"
#include "dyn_control_params.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{

CMATRIX compute_phase_corrections(CMATRIX& S, double tol);
CMATRIX compute_phase_corrections(CMATRIX& S);
vector<int> get_reordering(CMATRIX& time_overlap);
MATRIX make_cost_mat(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha);
vector<int> Munkres_Kuhn(CMATRIX& orb_mat_inp, CMATRIX& en_mat_inp, double alpha, int verbosity);
vector<int> get_stochastic_reordering(CMATRIX& time_overlap, Random& rnd);
vector<int> get_stochastic_reordering2(CMATRIX& time_overlap, Random& rnd);
vector<int> get_stochastic_reordering3(CMATRIX& time_overlap, Random& rnd, int convergence, int max_number_of_attempts);
vector<int> get_stochastic_reordering3(CMATRIX& time_overlap, Random& rnd, 
                                       int convergence, int max_number_of_attempts,
                                       double filter_tol, int verbosity_level
                                       );

vector<int> permute_states(vector<vector<int> >& perms, vector<int>& act_states);

CMATRIX permutation2cmatrix(vector<int>& permutation);
void update_projectors(dyn_control_params& prms, vector<CMATRIX>& projectors, 
  vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd);

vector< vector<int> > compute_permutations(dyn_control_params& prms, vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd);
vector<CMATRIX> compute_projectors(dyn_control_params& prms, vector<CMATRIX>& Eadi, vector<CMATRIX>& St, Random& rnd);
vector<CMATRIX> compute_projectors(dyn_control_params& prms, vector<CMATRIX>& St, vector<vector<int> >& perms);


CMATRIX raw_to_dynconsyst(CMATRIX& amplitudes, vector<CMATRIX>& projectors);
CMATRIX dynconsyst_to_raw(CMATRIX& amplitudes, vector<CMATRIX>& projectors);


}// namespace libdyn
}// liblibra

#endif // DYN_PROJECTORS_H

