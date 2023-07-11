/*********************************************************************************
* Copyright (C) 2015-2019 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file dyn_hop_acceptance.h
  \brief The header for dyn_hop_acceptance.cpp
    
*/

#ifndef DYN_HOP_ACCEPTANCE_H
#define DYN_HOP_ACCEPTANCE_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../nhamiltonian/libnhamiltonian.h"
#include "../io/libio.h"
#include "dyn_control_params.h"
#include "dyn_variables.h"
#include "../Units.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libnhamiltonian;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{

int can_rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, vector<int>& which_dofs);
int can_rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t);

void rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, int do_reverse, vector<int>& which_dofs);
void rescale_along_vector(double E_old, double E_new, MATRIX& p, MATRIX& invM, MATRIX& t, int do_reverse);

vector<double> Boltz_quant_prob(vector<double>& E, double T);
double Boltz_cl_prob(double E, double T);
double Boltz_cl_prob_up(double E, double T);
double HO_prob(vector<double>& E, vector<int>& qn, double T, vector<double>& prob);
double HO_prob_up(vector<double>& E, vector<int>& qn, double T, vector<double>& prob);
double boltz_factor(double E_new, double E_old, double T, int boltz_opt);


vector<int> accept_hops(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors, */
       nHamiltonian& ham, vector<int>& proposed_states, vector<int>& initial_states, Random& rnd,
       vector<int>& which_trajectories);

vector<int> accept_hops(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors, */
       nHamiltonian& ham, vector<int>& proposed_states, vector<int>& initial_states, Random& rnd);



vector<int> where_can_we_hop(int traj, dyn_control_params& prms,
       MATRIX& q, MATRIX& p,  MATRIX& invM, CMATRIX& Coeff, /*vector<CMATRIX>& projectors, */
       nHamiltonian& ham, vector<int>& act_states, Random& rnd);



void handle_hops_nuclear(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, /*vector<CMATRIX>& projectors,*/
       nHamiltonian& ham, vector<int>& new_states, vector<int>& old_states);



}// namespace libdyn
}// liblibra

#endif // DYN_HOP_ACCEPTANCE_H

