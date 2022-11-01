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


/// liblibra namespace
namespace liblibra{

using namespace libio;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{



///================  In dyn_methods_dish.cpp  ===================================

vector<int> decoherence_event(MATRIX& coherence_time, MATRIX& coherence_interval, int decoherence_event_option, Random& rnd);
vector<int> decoherence_event(MATRIX& coherence_time, MATRIX& coherence_interval, Random& rnd);

/*
vector<int> dish(dyn_control_params& prms,
       MATRIX& q, MATRIX& p,  MATRIX& invM, CMATRIX& Coeff, 
       nHamiltonian& ham, vector<int>& act_states, MATRIX& coherence_time, 
       vector<MATRIX>& decoherence_rates, Random& rnd);
*/

vector<int> dish_hop_proposal(vector<int>& act_states, CMATRIX& Coeff, 
  MATRIX& coherence_time, vector<MATRIX>& decoherence_rates, Random& rnd);

void dish_project_out_collapse(vector<int>& old_states, vector<int>& proposed_states, vector<int>& new_states, 
  CMATRIX& Coeff, MATRIX& coherence_time, int collapse_option);


/*
int ida(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi);
int dish(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib,
          int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2);
int dish(Electronic& el, Nuclear& mol, Hamiltonian& ham, 
          MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2);
*/




}// namespace libdyn
}// liblibra

#endif // DYN_METHODS_H

