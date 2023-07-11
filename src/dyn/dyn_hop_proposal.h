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
  \file dyn_hop_proposal.h
  \brief The header for dyn_hop_proposal.cpp
    
*/

#ifndef DYN_HOP_PROPOSAL_H
#define DYN_HOP_PROPOSAL_H

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

MATRIX hopping_probabilities_fssh(dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib);
vector<double> hopping_probabilities_fssh(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int act_state_indx);

MATRIX hopping_probabilities_gfsh(dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib);
vector<double> hopping_probabilities_gfsh(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int atc_state_indx);

vector<double> hopping_probabilities_fssh2(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& denmat_old, int act_state_indx);

MATRIX hopping_probabilities_mssh(dyn_control_params& prms, CMATRIX& Coeff, CMATRIX& Hvib);
vector<double> hopping_probabilities_mssh(dyn_control_params& prms, CMATRIX& denmat, CMATRIX& Hvib, int atc_state_indx);

//MATRIX compute_hopping_probabilities_lz(nHamiltonian* ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia);
//MATRIX compute_hopping_probabilities_lz(nHamiltonian& ham, int rep, MATRIX& p, const MATRIX& invM, MATRIX& prev_ham_dia);
vector<double> hopping_probabilities_lz(nHamiltonian* ham, nHamiltonian* ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM);
vector<double> hopping_probabilities_lz(nHamiltonian& ham, nHamiltonian& ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM);


vector<double> hopping_probabilities_zn(nHamiltonian* ham, nHamiltonian* ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM);
vector<double> hopping_probabilities_zn(nHamiltonian& ham, nHamiltonian& ham_prev, int act_state_indx, int rep, MATRIX& p, const MATRIX& invM);


vector<double> hopping_probabilities_mash(dyn_control_params& prms, CMATRIX& denmat);

/*
vector<MATRIX> hop_proposal_probabilities(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors,
       nHamiltonian& ham, vector<MATRIX>& prev_ham_dia);
*/
vector<MATRIX> hop_proposal_probabilities(dyn_control_params& prms,
       MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C,
       nHamiltonian& ham, vector<MATRIX>& prev_ham_dia);

vector< vector<double> > hop_proposal_probabilities(dyn_control_params& prms,
       dyn_variables& dyn_var, nHamiltonian& ham, nHamiltonian& ham_prev);

//vector< vector<double> > hop_proposal_probabilities(dyn_control_params& prms, 
//       dyn_variables& dyn_var, nHamiltonian& ham, vector<MATRIX>& prev_ham_dia);


int hop(vector<double>& prob, double ksi);
int hop(int initstate, MATRIX& g, double ksi);
int hop(int initstate, vector<double>& g, double ksi);

vector<int> propose_hops(vector<MATRIX>& g, vector<int>& act_states, Random& rnd);
vector<int> propose_hops(vector< vector<double> >& g, vector<int>& act_states, Random& rnd);




}// namespace libdyn
}// liblibra

#endif // DYN_HOP_PROPOSAL_H

