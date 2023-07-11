/*********************************************************************************
* Copyright (C) 2015-2021 Alexey V. Akimov
*
* This file is distributed under the terms of the GNU General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
* See the file LICENSE in the root directory of this distribution
* or <http://www.gnu.org/licenses/>.
*
*********************************************************************************/
/**
  \file Dynamics_Nuclear.h
  \brief The file describes the functions for nuclear (classical) dynamics
    
*/

#ifndef DYNAMICS_H
#define DYNAMICS_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../nhamiltonian/libnhamiltonian.h"
#include "../io/libio.h"
#include "thermostat/Thermostat.h"

#include "dyn_decoherence.h"
#include "dyn_control_params.h"
#include "dyn_hop_acceptance.h"
#include "dyn_hop_proposal.h"
#include "dyn_methods.h"
#include "dyn_projectors.h"
#include "dyn_variables.h"
#include "dyn_ham.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libnhamiltonian;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{

using namespace libthermostat;



//========== Dynamics.cpp ===================

void aux_get_transforms(CMATRIX** Uprev, nHamiltonian& ham);


// Adding the NBRA flag to the functions in the header
vector<CMATRIX> compute_St(nHamiltonian* ham, int isNBRA);
vector<CMATRIX> compute_St(nHamiltonian& ham, int isNBRA);
vector<CMATRIX> compute_St(nHamiltonian& ham);

vector<CMATRIX> compute_St(nHamiltonian* ham, nHamiltonian* ham_prev, int isNBRA);
vector<CMATRIX> compute_St(nHamiltonian& ham, nHamiltonian& ham_prev, int isNBRA);
vector<CMATRIX> compute_St(nHamiltonian& ham, nHamiltonian& ham_prev);


MATRIX momenta_on_excited_states(dyn_variables& dyn_var, nHamiltonian* ham, int itraj);
MATRIX momenta_on_excited_states(dyn_variables& dyn_var, nHamiltonian& ham, int itraj);

void SSY_correction(CMATRIX& Ham, dyn_variables& dyn_var, nHamiltonian* ham, int itraj);
void SSY_correction(CMATRIX& Ham, dyn_variables& dyn_var, nHamiltonian& ham, int itraj);

CMATRIX Zhu_Liouvillian(double Etot, CMATRIX& Ham, CMATRIX& rho);

void propagate_electronic(dyn_variables& dyn_var, nHamiltonian& ham, nHamiltonian& ham_prev, dyn_control_params& prms);
void propagate_electronic(dyn_variables& dyn_var, nHamiltonian* ham, nHamiltonian* ham_prev, dyn_control_params& prms);

/*
void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states, 
              nHamiltonian& ham, bp::object py_funct, bp::dict model_params, bp::dict dyn_params, Random& rnd);

void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states, 
              nHamiltonian& ham, bp::object py_funct, bp::dict& model_params, bp::dict& dyn_params, Random& rnd, 
              vector<Thermostat>& therm);

void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<CMATRIX>& projectors, vector<int>& act_states, 
              nHamiltonian& ham, bp::object py_funct, bp::dict& model_params, bp::dict& dyn_params, Random& rnd, 
              vector<Thermostat>& therm, dyn_variables& dyn_var);
*/

void compute_dynamics(dyn_variables& dyn_var, bp::dict dyn_params, nHamiltonian& ham, nHamiltonian& ham_aux, 
                      bp::object py_funct, bp::dict model_params, Random& rnd, vector<Thermostat>& therm);





}// namespace libdyn
}// liblibra

#endif // DYNAMICS_H

