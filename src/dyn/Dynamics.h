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
  \file Dynamics_Nuclear.h
  \brief The file describes the functions for nuclear (classical) dynamics
    
*/

#ifndef DYNAMICS_H
#define DYNAMICS_H

// External dependencies
#include "../math_linalg/liblinalg.h"
#include "../hamiltonian/libhamiltonian.h"
#include "../io/libio.h"
#include "thermostat/Thermostat.h"
#include "dyn_control_params.h"


/// liblibra namespace
namespace liblibra{

using namespace libio;
using namespace libhamiltonian;
namespace bp = boost::python;

/// libdyn namespace
namespace libdyn{

using namespace libthermostat;


// Verlet.cpp
void Verlet0(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params);
void Verlet1(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params);
void Verlet1(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, int entanglement_opt);

// Verlet_nvt.cpp
void Verlet0_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, Thermostat& therm);
void Verlet1_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, vector<Thermostat>& therm);
void Verlet1_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, int entanglement_opt, vector<Thermostat>& therm);
void Verlet1_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, Thermostat& therm);
void Verlet1_nvt(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, nHamiltonian& ham, bp::object py_funct, bp::object params, int entanglement_opt, Thermostat& therm);


// Ehrenfest.cpp
void Ehrenfest0(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, nHamiltonian& ham, bp::object py_funct, bp::object params, int rep);
void Ehrenfest1(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, nHamiltonian& ham, bp::object py_funct, bp::object params, int rep);
void Ehrenfest2(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, 
                nHamiltonian& ham, bp::object py_funct, bp::object params, int rep, int do_reordering, int do_phase_correction);
void Ehrenfest2(double dt, MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, 
                nHamiltonian& ham, bp::object py_funct, bp::object params, int rep);

//========== Dynamics.cpp ===================

MATRIX aux_get_forces(dyn_control_params& prms, nHamiltonian& ham, vector<int>& act_states, CMATRIX& amplitudes);
MATRIX aux_get_forces(bp::dict params, nHamiltonian& ham, vector<int>& act_states, CMATRIX& amplitudes);

void aux_get_transforms(CMATRIX** Uprev, nHamiltonian& ham);

void do_reordering(dyn_control_params& prms, nHamiltonian& ham,
                   vector<int>& act_states, CMATRIX& C, CMATRIX** Uprev, Random& rnd);


void do_phase_correction(dyn_control_params& prms, nHamiltonian& ham, 
                         vector<int>& act_states, CMATRIX& C, CMATRIX** Uprev);


void update_Hamiltonian_q(dyn_control_params& prms, MATRIX& q, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params);
void update_Hamiltonian_q(bp::dict prms, MATRIX& q, nHamiltonian& ham, 
                          bp::object py_funct, bp::object model_params);


void update_Hamiltonian_p(dyn_control_params& prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM);
void update_Hamiltonian_p(bp::dict prms, nHamiltonian& ham, MATRIX& p, MATRIX& invM);


CMATRIX transform_amplitudes(int rep_in, int rep_out, CMATRIX& C, nHamiltonian& ham);


void do_surface_hopping(dyn_control_params& prms,
              MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, vector<MATRIX>& prev_ham_dia, Random& rnd);
void do_surface_hopping(bp::dict prms,
              MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, vector<MATRIX>& prev_ham_dia, Random& rnd);



void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, bp::object py_funct, bp::dict model_params, bp::dict dyn_params, Random& rnd);



}// namespace libdyn
}// liblibra

#endif // DYNAMICS_H

