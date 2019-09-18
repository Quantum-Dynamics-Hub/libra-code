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

// Dynamics.coo
void compute_dynamics(MATRIX& q, MATRIX& p, MATRIX& invM, CMATRIX& C, vector<int>& act_states,
              nHamiltonian& ham, bp::object py_funct, bp::object params, boost::python::dict params1, Random& rnd);



}// namespace libdyn
}// liblibra

#endif // DYNAMICS_H

