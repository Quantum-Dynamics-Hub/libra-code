/*********************************************************************************
* Copyright (C) 2015-2017 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 2 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_decoherence.h
  \brief The file describes the functions for decoherence corrections
    
*/

#ifndef DYN_DECOHERENCE_H
#define DYN_DECOHERENCE_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../io/libio.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{


///================  In dyn_decoherence_methods.cpp  ===================================

CMATRIX sdm(CMATRIX& Coeff, double dt, int act_st, MATRIX& decoh_rates);
CMATRIX sdm(CMATRIX& Coeff, double dt, vector<int>& act_st, vector<MATRIX>& decoh_rates);

void project_out(CMATRIX& Coeff, int traj, int i);
void collapse(CMATRIX& Coeff, int traj, int i, int collapse_option);

void instantaneous_decoherence(CMATRIX& Coeff, 
   vector<int>& accepted_states, vector<int>& proposed_states, vector<int>& initial_states,
   int instantaneous_decoherence_variant, int collapse_option);



///================  In dyn_decoherence_time.cpp  ===================================

MATRIX edc_rates(CMATRIX& Hvib, double Ekin, double C_param, double eps_param);
vector<MATRIX> edc_rates(vector<CMATRIX>& Hvib, vector<double>& Ekin, double C_param, double eps_param);

void dephasing_informed_correction(MATRIX& decoh_rates, CMATRIX& Hvib, MATRIX& ave_gaps);
void dephasing_informed_correction(vector<MATRIX>& decoh_rates, vector<CMATRIX>& Hvib, MATRIX& ave_gaps);

MATRIX coherence_intervals(CMATRIX& Coeff, MATRIX& rates);


///================  In dyn_methods_dish.cpp  ===================================
/*
int ida(CMATRIX& Coeff, int old_st, int new_st, double E_old, double E_new, double T, double ksi);
int dish(Electronic& el, MATRIX& t_m, const MATRIX& tau_m, const CMATRIX& Hvib,
          int use_boltz_flag, double Ekin, double T, double ksi1, double ksi2);
int dish(Electronic& el, Nuclear& mol, Hamiltonian& ham, 
          MATRIX& t_m, const MATRIX& tau_m, int use_boltz_flag, double T, double ksi1, double ksi2);
*/




}// namespace libdyn
}// liblibra

#endif // DYN_DECOHERENCE_H

